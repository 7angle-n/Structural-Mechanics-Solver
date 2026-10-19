import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from fastapi import FastAPI, HTTPException, Request
from solver import solve, deflection, indeterminate, methods

app = FastAPI(title="Structural Mechanics Solver API")


@app.post("/api/solve")
@app.post("/api/index")
@app.post("/{path:path}")
async def solve_beam(req: Request):
    try:
        b = await req.json()
        if b.get("mode") == "deflection":
            return deflection(b, b.get("EI") or 1, b.get("at"))
        if b.get("mode") == "methods":
            return methods(b, b.get("EI") or 1, b.get("at"))
        if b.get("mode") == "indeterminate":
            return indeterminate(b, b.get("EI") or 1)
        return solve(b)
    except (ValueError, KeyError, TypeError) as e:
        raise HTTPException(400, f"Invalid problem: {e}")


if not os.environ.get("VERCEL"):  # local dev only: Vercel serves /public itself
    from fastapi.staticfiles import StaticFiles

    class NoCache(StaticFiles):  # always serve the newest files while developing
        async def get_response(self, path, scope):
            r = await super().get_response(path, scope)
            r.headers["Cache-Control"] = "no-store"
            return r

    app.mount("/", NoCache(directory=os.path.join(ROOT, "public"), html=True))
