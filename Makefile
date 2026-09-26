PY ?= python3

.PHONY: data site check clean

data:
	$(PY) code/make_pencil.py
	$(PY) code/solve_lines.py
	$(PY) code/line_data.py
	$(PY) code/export_data.py

site:
	$(PY) code/build_site.py

check: site
	$(PY) code/check_render.py

clean:
	rm -rf code/build
