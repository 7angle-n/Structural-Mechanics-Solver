# Structural Mechanics Solver

Created by Nafis Al Zahid.

Interactive web app: define a determinate beam (pin/roller or fixed supports; point loads, UDL, UVL, couples) and get reactions, section equations V(x) and M(x) with step-by-step working, numeric key values, and SFD/BMD plots.

**Stack:** static frontend (`public/`, KaTeX + Plotly via CDN) and a FastAPI function (`api/index.py`) using SymPy (`solver/core.py`).

## Run locally
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt
pytest                              # verify the solver
uvicorn api.index:app --reload      # open http://127.0.0.1:8000
```

## Deploy on Vercel
```bash
npm i -g vercel
vercel          # preview
vercel --prod
```
Or push to GitHub and import the repo in Vercel (framework preset: Other). `vercel.json` rewrites `/api/*` to the FastAPI app.

## Conventions
Applied loads down +, reactions up +, couples clockwise +, fixed-support moment M_A counter-clockwise +, sagging M +. The BMD is plotted with sagging downward (tension side).

## API
`POST /api/solve`
```json
{"length":10,
 "supports":[{"type":"pin","x":0},{"type":"roller","x":10}],
 "loads":[{"type":"point","x":3,"magnitude":10},
          {"type":"udl","x_start":5,"x_end":10,"w":4},
          {"type":"uvl","x_start":0,"x_end":4,"w_start":0,"w_end":6},
          {"type":"moment","x":6,"magnitude":20}]}
```

## Pages
| Path | What it does | Runs on |
|---|---|---|
| `/` | Dark engineering-lab homepage: live 3D truss (direct stiffness solve), playground, capability map, results workspace, topic wheel | static |
| `/docs/` | Conventions, run and deploy notes | static |
| `/beam/` | SFD/BMD, reactions; deflection by double integration, conjugate beam, moment-area, superposition and unit load; indeterminate beams (compatibility) | `/api/solve` (SymPy) |
| `/inertia/` | Centroids, I, parallel axis, composite sections | browser |
| `/mohr/` | Mohr's circle, principal stresses | browser |
| `/axial/`, `/torsion/` | Stepped bars (thermal strain, gaps, fixed-fixed compatibility); shafts, gear pairs, shaft sizing from power | browser |
| `/truss/` | Determinate trusses: method of joints and method of sections | browser |
| `/columns/` | Euler, Johnson, secant | browser |
| `/stress/` | Bending and shear stress distribution | browser |

`POST /api/solve` accepts `"mode": "deflection" | "methods" | "indeterminate"` (default: determinate SFD/BMD).

Open the site through `http://127.0.0.1:8000` (or Vercel). All page paths are relative, so the static pages also work from any static server rooted at `public/`; only `/beam/` needs the API.

## Fonts
Text uses Latin Modern (the LaTeX font) and Computer Modern Typewriter, bundled in `public/assets/fonts` (no Google Fonts or CDN needed for fonts). Equations are rendered by KaTeX.
