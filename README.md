# Mr. Imputation

Mr. Imputation is an educational synthetic genotype-imputation explorer.
The frontend is React/TypeScript/Vite, and the statistical backend is
Python/FastAPI.

## Local development

Start the backend from the project root:

```bash
PYTHONPATH=. python3 -m uvicorn backend.app.api:app --reload
```

In another terminal, start the frontend:

```bash
cd frontend
npm install
npm run dev
```

The Vite development proxy sends `/api` requests to
`http://127.0.0.1:8000`. The browser application is available at
`http://localhost:5173`.

## Free public deployment

The application uses two services:

1. A Render free web service for FastAPI.
2. A Cloudflare Pages static site for the built React frontend.

### Backend on Render

Create a web service from this repository, using:

- Build command: `pip install -r backend/requirements.txt`
- Start command: `uvicorn backend.app.api:app --host 0.0.0.0 --port $PORT`
- Health check path: `/health`

Set `FRONTEND_ORIGIN` to the final Cloudflare Pages URL, for example
`https://mr-imputation.pages.dev`.

The included `render.yaml` contains the service definition. Render free
services can sleep when idle, so the first request after inactivity may be
slower.

### Frontend on Cloudflare Pages

Create a Pages project from the same repository with:

- Root directory: `frontend`
- Build command: `npm run build`
- Output directory: `dist`
- Environment variable:
  - `VITE_API_BASE_URL=https://<your-render-service>.onrender.com`

The frontend falls back to the local `/api` Vite proxy when this variable is
unset. Do not put credentials in frontend environment variables.

After both services are deployed, open the Cloudflare Pages URL. The browser
will call the Render API rather than a `127.0.0.1` address.
