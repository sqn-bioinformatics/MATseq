"""Lightweight disk-based caching system for pipeline steps and GO files."""

import pickle
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Optional, Dict
from collections import namedtuple


class PipelineCache:
    """Manage caching of intermediate results for reproducibility and speed."""

    def __init__(self, cache_dir: Path = None):
        """Initialize cache manager.

        Args:
            cache_dir: Directory for storing cached files. Defaults to results/cache.
        """
        if cache_dir is None:
            cache_dir = Path.cwd() / "results" / "cache"

        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.cache_dir / "manifest.json"
        self.manifest = self._load_manifest()

    def _load_manifest(self) -> dict:
        """Load manifest of cached files.

        Returns:
            Dictionary mapping cache keys to file metadata.
        """
        if self.manifest_path.exists():
            try:
                with open(self.manifest_path, "r") as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Warning: Could not load manifest, starting fresh: {e}")
                return {}
        return {}

    def _save_manifest(self):
        """Save manifest of cached files."""
        try:
            with open(self.manifest_path, "w") as f:
                json.dump(self.manifest, f, indent=2)
        except IOError as e:
            print(f"Warning: Could not save manifest: {e}")

    def _get_cache_key(self, name: str, params: dict = None) -> str:
        """Generate cache key from name and parameters.

        Args:
            name: Name of the cached item.
            params: Dictionary of parameters affecting the cache.

        Returns:
            Hash-based cache key.
        """
        key_str = name
        if params:
            # Create deterministic hash of parameters
            params_json = json.dumps(params, sort_keys=True)
            param_hash = hashlib.md5(params_json.encode()).hexdigest()[:8]
            key_str = f"{name}_{param_hash}"

        return key_str

    def _get_cache_path(self, cache_key: str) -> Path:
        """Get file path for cache key.

        Args:
            cache_key: Cache key identifier.

        Returns:
            Path to cached file.
        """
        return self.cache_dir / f"{cache_key}.pkl"

    def get(self, name: str, params: dict = None) -> Optional[Any]:
        """Retrieve cached item if it exists.

        Args:
            name: Name of the cached item.
            params: Parameters used to generate the item (for cache validity).

        Returns:
            Cached object if available, None otherwise.
        """
        cache_key = self._get_cache_key(name, params)
        cache_path = self._get_cache_path(cache_key)

        if cache_path.exists() and cache_key in self.manifest:
            try:
                with open(cache_path, "rb") as f:
                    return pickle.load(f)
            except Exception as e:
                print(f"Warning: Could not load cache {cache_key}: {e}")
                return None

        return None

    def set(
        self, name: str, obj: Any, params: dict = None, description: str = ""
    ) -> Path:
        """Cache an object to disk.

        Args:
            name: Name of the cached item.
            obj: Object to cache (must be picklable).
            params: Parameters used to generate the item.
            description: Human-readable description of cached content.

        Returns:
            Path to cached file.
        """
        cache_key = self._get_cache_key(name, params)
        cache_path = self._get_cache_path(cache_key)

        try:
            with open(cache_path, "wb") as f:
                pickle.dump(obj, f)

            self.manifest[cache_key] = {
                "name": name,
                "description": description,
                "params": params,
                "path": str(cache_path),
            }
            self._save_manifest()

            return cache_path
        except Exception as e:
            print(f"Warning: Could not cache {name}: {e}")
            return cache_path

    def clear(self, name: str = None, params: dict = None):
        """Clear cache entries.

        Args:
            name: Specific cache name to clear. If None, clears all.
            params: Parameters for partial key matching.
        """
        if name is None:
            # Clear all
            for cache_key in list(self.manifest.keys()):
                cache_path = self._get_cache_path(cache_key)
                if cache_path.exists():
                    cache_path.unlink()
            self.manifest = {}
        else:
            cache_key = self._get_cache_key(name, params)
            cache_path = self._get_cache_path(cache_key)
            if cache_path.exists():
                cache_path.unlink()
            if cache_key in self.manifest:
                del self.manifest[cache_key]

        self._save_manifest()

    def cached_call(
        self,
        func: Callable,
        name: str,
        params: dict = None,
        *args,
        force_recompute: bool = False,
        **kwargs,
    ) -> Any:
        """Call a function with caching.

        Args:
            func: Function to call.
            name: Cache name for this computation.
            params: Parameters dictionary for cache key.
            force_recompute: Skip cache and recompute.
            *args: Positional arguments for func.
            **kwargs: Keyword arguments for func.

        Returns:
            Result of function call (from cache or fresh computation).
        """
        if not force_recompute:
            cached_result = self.get(name, params)
            if cached_result is not None:
                print(f"Loading {name} from cache...")
                return cached_result

        print(f"Computing {name}...")
        result = func(*args, **kwargs)
        self.set(name, result, params, description=f"Result of {func.__name__}")

        return result

    def list_cached(self) -> list:
        """List all cached items.

        Returns:
            List of cached item metadata.
        """
        return list(self.manifest.values())


# Global memory cache for GO files
_go_dag_cache: Optional[Any] = None
_gene2go_cache: Optional[Any] = None
_geneid2nt_cache: Optional[Dict[int, Any]] = None


def load_geneid2nt(gene_file: Path) -> Dict[int, Any]:
    """Load GENEID2NT from gene_result_ncbi_human_proteincoding.txt file (cached).

    Args:
        gene_file: Path to gene file.

    Returns:
        Dictionary mapping GeneID to ntncbi namedtuples.
    """
    global _geneid2nt_cache
    if _geneid2nt_cache is not None:
        return _geneid2nt_cache

    ntncbi = namedtuple(
        "ntncbi",
        "tax_id Org_name GeneID CurrentID Status Symbol Aliases description other_designations map_location chromosome genomic_nucleotide_accession_version start_position_on_the_genomic_accession end_position_on_the_genomic_accession orientation exon_count OMIM",
    )

    geneid2nt = {}
    try:
        with open(gene_file, "r", encoding="utf-8") as f:
            next(f)
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) >= 17:
                    try:
                        gene_id = int(fields[2])
                        genomic_acc = fields[11].replace(".", "_")
                        entry = ntncbi(
                            tax_id=fields[0],
                            Org_name=fields[1],
                            GeneID=gene_id,
                            CurrentID=fields[3],
                            Status=fields[4],
                            Symbol=fields[5],
                            Aliases=fields[6],
                            description=fields[7],
                            other_designations=fields[8],
                            map_location=fields[9],
                            chromosome=fields[10],
                            genomic_nucleotide_accession_version=genomic_acc,
                            start_position_on_the_genomic_accession=fields[12],
                            end_position_on_the_genomic_accession=fields[13],
                            orientation=fields[14],
                            exon_count=fields[15],
                            OMIM=fields[16],
                        )
                        geneid2nt[gene_id] = entry
                    except (ValueError, IndexError):
                        continue
        _geneid2nt_cache = geneid2nt
        return geneid2nt
    except Exception as e:
        print(f"Warning: Could not load gene file {gene_file}: {e}")
        _geneid2nt_cache = {}
        return {}


def load_go_dag(obo_file: Path) -> Any:
    """Load GO DAG from obo file (cached).

    Args:
        obo_file: Path to go-basic.obo file.

    Returns:
        GODag instance.
    """
    global _go_dag_cache
    if _go_dag_cache is not None:
        return _go_dag_cache

    try:
        from goatools.obo_parser import GODag as _GODag
        _go_dag_cache = _GODag(str(obo_file))
        return _go_dag_cache
    except Exception as e:
        print(f"Warning: Could not load GO DAG from {obo_file}: {e}")
        return None


def load_gene2go(gene2go_file: Path, taxids: list = None, namespaces: set = None) -> Any:
    """Load gene to GO associations (cached).

    Args:
        gene2go_file: Path to gene2go file.
        taxids: List of taxonomy IDs to load (default: [9606] for human).
        namespaces: Set of GO namespaces (default: {'BP'}).

    Returns:
        Gene2GoReader instance.
    """
    global _gene2go_cache
    if _gene2go_cache is not None:
        return _gene2go_cache

    if taxids is None:
        taxids = [9606]
    else:
        taxids = [int(t) for t in taxids]
    if namespaces is None:
        namespaces = {"BP"}
    else:
        namespaces = set(namespaces)

    try:
        from goatools.anno.genetogo_reader import Gene2GoReader as _Gene2GoReader
        _gene2go_cache = _Gene2GoReader(
            str(gene2go_file), taxids=taxids, namespaces=namespaces
        )
        return _gene2go_cache
    except Exception as e:
        print(f"Warning: Could not load gene2go from {gene2go_file}: {e}")
        return None


def clear_go_cache() -> None:
    """Clear all cached GO data."""
    global _go_dag_cache, _gene2go_cache, _geneid2nt_cache
    _go_dag_cache = None
    _gene2go_cache = None
    _geneid2nt_cache = None
