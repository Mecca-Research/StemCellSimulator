// Minimal static file server for tests (no dependencies).
import http from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import { extname, join, normalize, resolve } from 'node:path';

const TYPES = {
  '.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8', '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8', '.json': 'application/json', '.png': 'image/png', '.jpg': 'image/jpeg',
  '.bin': 'application/octet-stream', '.svg': 'image/svg+xml', '.webp': 'image/webp', '.gif': 'image/gif',
};

export function serve(root, port = 0) {
  const base = resolve(root);
  const server = http.createServer(async (req, res) => {
    try {
      let path = decodeURIComponent(new URL(req.url, 'http://x').pathname);
      if (path.endsWith('/')) path += 'index.html';
      const file = normalize(join(base, path));
      if (!file.startsWith(base)) { res.writeHead(403).end(); return; }
      await stat(file);
      const body = await readFile(file);
      res.writeHead(200, { 'content-type': TYPES[extname(file)] ?? 'application/octet-stream', 'cache-control': 'no-store' });
      res.end(body);
    } catch {
      res.writeHead(404).end('not found');
    }
  });
  return new Promise((ok) => server.listen(port, '127.0.0.1', () => ok({ server, url: `http://127.0.0.1:${server.address().port}/` })));
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const { url } = await serve(process.argv[2] ?? '.', Number(process.argv[3] ?? 8080));
  console.log(`serving on ${url}`);
}
