/// <reference types="vitest/config" />
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// `npm run build` writes the app into the Python package, where `cellfactory-api` serves it at `/`.
// Fixed file names (no content hashes) keep the committed build's diffs small and reviewable.
// `npm run dev` serves on :5173 and forwards /api to a `cellfactory-api` on :8000.
export default defineConfig({
  plugins: [react()],
  base: "./",
  build: {
    outDir: "../src/cellfactory/web/static",
    emptyOutDir: true,
    sourcemap: false,
    rollupOptions: {
      output: {
        entryFileNames: "assets/app.js",
        chunkFileNames: "assets/[name].js",
        assetFileNames: "assets/[name][extname]",
      },
    },
  },
  server: {
    proxy: { "/api": "http://127.0.0.1:8000" },
  },
  test: { environment: "node" },
});
