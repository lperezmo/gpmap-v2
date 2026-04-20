# gpmap-v2 streamlit demo

Interactive, page-per-topic tour of the `gpmap-v2` API. Half show-and-tell, half
documentation. Hosted at [gpmap-v2.streamlit.app](https://gpmap-v2.streamlit.app).

## Run locally

```bash
pip install -r examples/streamlit/requirements.txt
streamlit run examples/streamlit/showcase.py
```

From a local editable checkout:

```bash
uv sync
uv run maturin develop --release
uv pip install streamlit plotly
uv run streamlit run examples/streamlit/showcase.py
```

## Layout

```
examples/streamlit/
  showcase.py              entry point, navigation, header, footer
  app_pages/               one file per topic
  utils/                   shared widget + chart helpers
  requirements.txt         streamlit cloud dependencies
```

Structure mirrors the multi-page pattern used by
[streamlit-aggrid-v2](https://github.com/lperezmo/st-aggrid) (the reference app
this demo is based on).
