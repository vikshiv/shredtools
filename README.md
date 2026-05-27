# **shredtools browser**: extract syntenic regions from HPRC



This branch hosts a small, static, in-browser UI for running `shredtools extract` queries against precomputed multi-MUM indexes for pangenomes from the Human Pangenome Reference Consortium (HPRC). It runs entirely client-side via Pyodide and reads the index/MUM data from public S3 using HTTP Range requests.

**Browser link:** [https://vikshiv.github.io/shredtools/](https://vikshiv.github.io/shredtools/)

---

## Run locally

Serve this directory over HTTP:

```bash
python -m http.server 8000
```

Then open [http://localhost:8000](http://localhost:8000) in your browser.