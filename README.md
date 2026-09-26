# Lines in a quartic pencil

An interactive page showing a real pencil of quartic surfaces X_λ = {f − λQ = 0} in ℝℙ³.
You can sweep λ, stop on each real line of the pencil, rotate a plane H through the line,
and watch the nodes of the plane section H ∩ X_λ add up to the local index of the line,
marked at Q:

    ind_Q(L) = ⟨ −Res(F₁, F₂) · Σ_p Q(p) / (F_⊥(p) · F_H′(p)) ⟩

## Layout

    code/            computation and site build
      make_pencil.py   step 1: choose f and Q (seeded)
      solve_lines.py   step 2: all 320 lines by homotopy continuation
      line_data.py     step 3: per-line data and index signs, with consistency checks
      export_data.py   step 4: write code/data/pencil.json
      build_site.py    inline the data into template.html -> site/index.html
      template.html    the page (WebGL surface renderer, controls, inset, node table)
      check_render.py  optional headless-Chromium screenshot check
      data/pencil.json the precomputed data used by the site (committed)
      build/           intermediate files (ignored)
    site/            the built static site (index.html only)

## Building

The site only needs the committed data, so a deploy step is just

    python3 code/build_site.py

and then publish the `site/` folder to any static host (GitHub Pages, Netlify, ...).
It uses only the Python standard library.

To recompute the pencil and its lines from scratch (a few minutes; results are seeded and reproducible):

    pip install -r code/requirements.txt
    make data     # steps 1-4
    make site

## Notes

- The page is a single self-contained `index.html`. It loads fonts from Google Fonts and falls back to system fonts without them. WebGL is required for the 3D view.
- Signs use the Kass–Wickelgren chart S = span(e₀ + a e₂ + b e₃, e₁ + c e₂ + d e₃),
  the monomial basis u⁴, …, r⁴ and det[Q|_L, uF₁, uF₂, rF₁, rF₂]. They are checked for invariance
  under random changes of basis and between charts; a different global convention would flip all signs.
- For this pencil: 320 complex lines, 12 real, and the local indices over ℝ sum to −8.
