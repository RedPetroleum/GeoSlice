#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path

TOKEN_RE = re.compile(r'([A-Za-z])\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')


def parse_tokens(line: str):
    return [(m.group(1).upper(), float(m.group(2)), m.start(), m.end()) for m in TOKEN_RE.finditer(line)]


def rebuild_line(original: str, replacements: dict[str, float], decimals: int) -> str:
    out = []
    last = 0
    for m in TOKEN_RE.finditer(original):
        axis = m.group(1).upper()
        out.append(original[last:m.start()])
        if axis in replacements:
            out.append(f"{axis}{replacements[axis]:.{decimals}f}")
        else:
            out.append(original[m.start():m.end()])
        last = m.end()
    out.append(original[last:])
    return ''.join(out)


def convert_non_tcp_to_tcp(xf: float, yf: float, zf: float, uf: float, vf: float, h: float, model: str, r: float):
    # Non-TCP output in this codebase does:
    # 1) swap X/Y before writing
    # 2) negate B/C before writing
    # Therefore recover machine-space first:
    mx = yf
    my = xf
    mz = zf

    # Recover internal B/C (deg)
    b_deg = -uf
    c_deg = -vf

    b = math.radians(b_deg)
    c = math.radians(c_deg)

    if model == "old":
        # old _getXYZ (rotationary tilt table)
        # mx = px*cos(c) - py*sin(c)
        # my = px*sin(c) + py*cos(c) - h*sin(b)
        # mz = pz - h*(1-cos(b))
        s = my + h * math.sin(b)
        px = mx * math.cos(c) + s * math.sin(c)
        py = -mx * math.sin(c) + s * math.cos(c)
        pz = mz + h * (1.0 - math.cos(b))
    else:
        # newconfig _getXYZ_newConfig (rotationary tilt table)
        # mx = px - r*cos(c) + h*sin(c)*sin(b) + r
        # my = py - r*sin(c) - h*cos(c)*sin(b)
        # mz = pz + h*cos(b) - h
        px = mx + r * math.cos(c) - h * math.sin(c) * math.sin(b) - r
        py = my + r * math.sin(c) + h * math.cos(c) * math.sin(b)
        pz = mz - h * math.cos(b) + h

    # TCP gcode in this project writes internal values directly (no swap/sign flip)
    return px, py, pz, b_deg, c_deg


def main():
    p = argparse.ArgumentParser(description="Convert non-TCP GCode (from this project) back to TCP-style GCode.")
    p.add_argument("input", type=Path, help="Input non-TCP gcode file")
    p.add_argument("output", type=Path, help="Output TCP gcode file")
    p.add_argument("--h", type=float, required=True, help="Tool length h used during generation")
    p.add_argument("--model", choices=["old", "newconfig"], default="old", help="Kinematic model used in generation")
    p.add_argument("--r", type=float, default=0.0, help="Radius r (required for --model newconfig)")
    p.add_argument("--decimals", type=int, default=3, help="Decimal places for rewritten axes")
    args = p.parse_args()

    if args.model == "newconfig" and abs(args.r) < 1e-12:
        print("Warning: --model newconfig uses parameter r; currently r is 0.0")

    lines = args.input.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    out_lines = []

    for line in lines:
        s = line.strip()
        if not s or s.startswith(";") or s.startswith("("):
            out_lines.append(line)
            continue

        tokens = parse_tokens(line)
        vals = {a: v for a, v, _, _ in tokens}

        # Convert only motion lines carrying full 5-axis pose
        if ("G" in vals) and int(vals["G"]) == 1 and all(k in vals for k in ("X", "Y", "Z", "U", "V")):
            px, py, pz, b_deg, c_deg = convert_non_tcp_to_tcp(
                vals["X"], vals["Y"], vals["Z"], vals["U"], vals["V"], args.h, args.model, args.r
            )
            repl = {"X": px, "Y": py, "Z": pz, "U": b_deg, "V": c_deg}
            out_lines.append(rebuild_line(line, repl, args.decimals))
        else:
            out_lines.append(line)

    args.output.write_text(''.join(out_lines), encoding="utf-8")
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
