#!/usr/bin/env python3
"""
topo.py - Get the topography of any latitude/longitude.

Samples an elevation grid around the point (via the free Open-Topo-Data API,
no API key needed) and reports:
  - elevation (m / ft)
  - slope (degrees and %), using Horn's method (same as GIS tools)
  - aspect (which way the ground faces, as a compass direction)
  - local relief and landform type (peak, ridge, slope, valley, flat...)
  - terrain steepness class
  - an ASCII elevation map of the surrounding area

Standard library only. Python 3.8+.

Usage:
  python topo.py 14.5995 120.9842
  python topo.py 27.9881 86.9250 --radius 1000 --step 100
  python topo.py 46.8523 -121.7603 --dataset ned10m      # US only, finer data
  python topo.py 14.5995 120.9842 --json
  python topo.py --selftest                              # offline math check
"""
import argparse
import json
import math
import sys
import time
import urllib.error
import urllib.request

API = "https://api.opentopodata.org/v1/{dataset}"
M_PER_DEG_LAT = 111_320.0


# ---------------------------------------------------------------- sampling
def build_grid(lat, lon, step, k):
    """Return list of rows (north -> south) of (lat, lon) tuples, (2k+1)^2 points."""
    coslat = max(math.cos(math.radians(lat)), 1e-6)
    rows = []
    for i in range(2 * k + 1):
        dy = (k - i) * step          # metres north of centre
        row = []
        for j in range(2 * k + 1):
            dx = (j - k) * step      # metres east of centre
            row.append((lat + dy / M_PER_DEG_LAT,
                        lon + dx / (M_PER_DEG_LAT * coslat)))
        rows.append(row)
    return rows


def fetch_elevations(points, dataset="srtm30m", retries=3):
    """Query Open-Topo-Data. points: list of (lat, lon). Returns list of floats."""
    out = []
    for s in range(0, len(points), 100):          # API limit: 100 per request
        chunk = points[s:s + 100]
        loc = "|".join(f"{la:.6f},{lo:.6f}" for la, lo in chunk)
        url = API.format(dataset=dataset) + "?locations=" + loc
        for attempt in range(retries):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "topo.py/1.0"})
                with urllib.request.urlopen(req, timeout=30) as r:
                    data = json.load(r)
                break
            except (urllib.error.URLError, TimeoutError) as e:
                if attempt == retries - 1:
                    raise SystemExit(f"Elevation API request failed: {e}")
                time.sleep(1.5 * (attempt + 1))
        if data.get("status") != "OK":
            raise SystemExit(f"API error: {data.get('error', data)}")
        for res in data["results"]:
            if res["elevation"] is None:
                raise SystemExit("No elevation data at this location for "
                                 f"dataset '{dataset}' (try --dataset aster30m "
                                 "or --dataset etopo1 for oceans/poles).")
            out.append(float(res["elevation"]))
        if s + 100 < len(points):
            time.sleep(1.1)                       # free tier: ~1 request/second
    return out


# ---------------------------------------------------------------- analysis
COMPASS = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]


def compass(deg):
    return COMPASS[int((deg + 22.5) % 360 // 45)]


def slope_aspect(z3, step):
    """Horn's method on a 3x3 window (rows north->south). Returns (slope_deg, aspect_deg|None)."""
    a, b, c = z3[0]
    d, _, f = z3[1]
    g, h, i = z3[2]
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * step)   # +ve = rises to east
    dzdy = ((a + 2 * b + c) - (g + 2 * h + i)) / (8 * step)   # +ve = rises to north
    slope = math.degrees(math.atan(math.hypot(dzdx, dzdy)))
    if math.hypot(dzdx, dzdy) < 1e-9:
        return slope, None
    aspect = math.degrees(math.atan2(-dzdx, -dzdy)) % 360     # direction of downslope
    return slope, aspect


def steepness_class(slope):
    for limit, name in [(2, "flat"), (5, "gently sloping"), (10, "moderately sloping"),
                        (15, "strongly sloping"), (30, "steep"), (45, "very steep")]:
        if slope < limit:
            return name
    return "extremely steep / cliff"


def landform(grid, slope):
    flat = [v for row in grid for v in row]
    k = len(grid) // 2
    centre = grid[k][k]
    relief = max(flat) - min(flat)
    others = [v for v in flat]
    tpi = centre - (sum(others) - centre) / (len(others) - 1)   # topographic position
    if relief < 3:
        return "flat / plain", relief, tpi
    r = tpi / relief
    # Slope is ~0 at a true summit or basin floor, so don't call it flat on slope alone.
    if slope < 1 and abs(r) <= 0.10:
        return "flat / plain", relief, tpi
    if r > 0.30:
        return ("peak / summit" if slope < 15 else "ridge / crest"), relief, tpi
    if r < -0.30:
        return ("depression / valley floor" if slope < 15 else "valley / gully"), relief, tpi
    if r > 0.10:
        return "upper slope / shoulder", relief, tpi
    if r < -0.10:
        return "lower slope / footslope", relief, tpi
    return "mid-slope", relief, tpi


def ascii_map(grid):
    shades = " .:-=+*#%@"
    flat = [v for row in grid for v in row]
    lo, hi = min(flat), max(flat)
    span = (hi - lo) or 1
    lines = []
    for row in grid:
        lines.append(" ".join(shades[int((v - lo) / span * (len(shades) - 1))] * 2 for v in row))
    return "\n".join(lines) + f"\n(N is up; ' '={lo:.0f} m ... '@'={hi:.0f} m)"


def analyze(lat, lon, grid, step):
    k = len(grid) // 2
    z3 = [grid[k - 1 + r][k - 1:k + 2] for r in range(3)]
    slope, aspect = slope_aspect(z3, step)
    form, relief, tpi = landform(grid, slope)
    z = grid[k][k]
    return {
        "latitude": lat, "longitude": lon,
        "elevation_m": round(z, 1), "elevation_ft": round(z * 3.28084, 1),
        "slope_deg": round(slope, 1), "slope_percent": round(math.tan(math.radians(slope)) * 100, 1),
        "aspect_deg": None if aspect is None else round(aspect, 0),
        "faces": "none (flat)" if aspect is None else compass(aspect),
        "steepness": steepness_class(slope),
        "landform": form,
        "local_relief_m": round(relief, 1),
        "area_min_m": round(min(v for r in grid for v in r), 1),
        "area_max_m": round(max(v for r in grid for v in r), 1),
        "sample_spacing_m": step,
        "area_size_m": step * (len(grid) - 1),
    }


# ---------------------------------------------------------------- main
def topography(lat, lon, radius=300, step=30, dataset="srtm30m", fetch=fetch_elevations):
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        raise SystemExit("Latitude must be -90..90 and longitude -180..180.")
    k = max(1, int(round(radius / step)))
    pts = build_grid(lat, lon, step, k)
    flat = [p for row in pts for p in row]
    elev = fetch(flat, dataset)
    n = 2 * k + 1
    grid = [elev[r * n:(r + 1) * n] for r in range(n)]
    return analyze(lat, lon, grid, step), grid


def selftest():
    # Plane rising to the east at 10% -> slope 5.71 deg, ground faces WEST (270).
    def fake(points, dataset):
        lat0, lon0 = points[len(points) // 2]
        coslat = math.cos(math.radians(lat0))
        return [100 + 0.10 * (lo - lon0) * M_PER_DEG_LAT * coslat for _, lo in points]
    res, _ = topography(14.6, 120.98, 150, 30, fetch=fake)
    assert abs(res["slope_deg"] - 5.7) < 0.2, res
    assert res["faces"] == "W", res
    # Plane rising to the north -> faces S.
    def fake2(points, dataset):
        lat0, _ = points[len(points) // 2]
        return [100 + 0.2 * (la - lat0) * M_PER_DEG_LAT for la, _ in points]
    res2, _ = topography(14.6, 120.98, 150, 30, fetch=fake2)
    assert res2["faces"] == "S" and abs(res2["slope_deg"] - 11.3) < 0.3, res2
    # Cone peak -> summit/ridge.
    def fake3(points, dataset):
        la0, lo0 = points[len(points) // 2]
        cl = math.cos(math.radians(la0))
        return [500 - 0.05 * math.hypot((la - la0) * M_PER_DEG_LAT, (lo - lo0) * M_PER_DEG_LAT * cl)
                for la, lo in points]
    res3, _ = topography(14.6, 120.98, 150, 30, fetch=fake3)
    assert "peak" in res3["landform"] or "ridge" in res3["landform"], res3
    print("selftest OK")


def main():
    ap = argparse.ArgumentParser(description="Topography of a latitude/longitude.")
    ap.add_argument("lat", nargs="?", type=float)
    ap.add_argument("lon", nargs="?", type=float)
    ap.add_argument("--radius", type=float, default=300, help="analysis radius in metres (default 300)")
    ap.add_argument("--step", type=float, default=30, help="grid spacing in metres (default 30; match dataset resolution)")
    ap.add_argument("--dataset", default="srtm30m",
                    help="srtm30m (global 60N-56S), aster30m (global), ned10m (US), mapzen, etopo1 (coarse, incl. ocean)")
    ap.add_argument("--json", action="store_true", help="output JSON")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        return selftest()
    if a.lat is None or a.lon is None:
        ap.error("lat and lon are required")

    res, grid = topography(a.lat, a.lon, a.radius, a.step, a.dataset)
    if a.json:
        print(json.dumps(res, indent=2))
        return
    print(f"Location: {res['latitude']}, {res['longitude']}  (dataset: {a.dataset})")
    print(f"Elevation : {res['elevation_m']} m  ({res['elevation_ft']} ft)")
    print(f"Slope     : {res['slope_deg']}°  ({res['slope_percent']}%)  - {res['steepness']}")
    print(f"Aspect    : faces {res['faces']}" + ("" if res['aspect_deg'] is None else f" ({res['aspect_deg']:.0f}°)"))
    print(f"Landform  : {res['landform']}")
    print(f"Relief    : {res['local_relief_m']} m within {res['area_size_m']:.0f} m x {res['area_size_m']:.0f} m "
          f"({res['area_min_m']}–{res['area_max_m']} m)")
    print("\nElevation map:")
    print(ascii_map(grid))


if __name__ == "__main__":
    main()