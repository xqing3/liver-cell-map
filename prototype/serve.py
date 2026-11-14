"""Serve the prototype with gzip, which cuts the full data from about 9.4 MB to 3.8 MB on the wire.

Usage (from the prototype folder):  python3 serve.py [port]     (default 8000)

python3 -m http.server also works but sends files uncompressed, so the page loads slower on a
real network. Compressed copies are kept in memory and rebuilt when a file changes.
"""
import gzip
import sys
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

COMPRESS = {".html", ".json", ".txt", ".js", ".css", ".csv"}
_cache = {}


class GzipHandler(SimpleHTTPRequestHandler):
    def send_head(self):
        path = Path(self.translate_path(self.path))
        if path.is_dir():
            path = path / "index.html"
        if path.suffix not in COMPRESS or not path.is_file() or "gzip" not in self.headers.get("Accept-Encoding", ""):
            return super().send_head()
        stamp = path.stat().st_mtime_ns
        if _cache.get(path, (None,))[0] != stamp:
            _cache[path] = (stamp, gzip.compress(path.read_bytes(), 6))
        body = _cache[path][1]
        self.send_response(200)
        self.send_header("Content-Type", self.guess_type(str(path)))
        self.send_header("Content-Encoding", "gzip")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Vary", "Accept-Encoding")
        self.end_headers()
        if self.command == "GET":
            self.wfile.write(body)
        return None


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    handler = partial(GzipHandler, directory=str(Path(__file__).parent))
    print(f"Serving the prototype at http://localhost:{port}/ (gzip on)")
    ThreadingHTTPServer(("", port), handler).serve_forever()
