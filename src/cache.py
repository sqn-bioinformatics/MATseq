"""Lightweight disk-based caching system for pipeline steps."""

import pickle
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Optional


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
            with open(self.manifest_path, "r") as f:
                return json.load(f)
        return {}

    def _save_manifest(self):
        """Save manifest of cached files."""
        with open(self.manifest_path, "w") as f:
            json.dump(self.manifest, f, indent=2)

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
