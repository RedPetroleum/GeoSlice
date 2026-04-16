import re
import math
import sys

def parse_gcode_line(line):
    if not line.startswith("G1"):
        return None
    vals = {}
    for token in line.split():
        if len(token) > 1 and token[0] in "XYZUVEF":
            try:
                vals[token[0]] = float(token[1:])
            except ValueError:
                pass
    if not all(k in vals for k in ("X", "Y", "Z", "U", "V")):
        return None
    return vals

def to_print_position(mx, my, mz, rad_b, rad_c):
    sinB, cosB = math.sin(rad_b), math.cos(rad_b)
    sinC, cosC = math.sin(rad_c), math.cos(rad_c)

    # Undo the X/Y swap from ToMachinePosition
    X_neu = my
    Y_neu = mx
    Z_neu = mz

    # Apply M^T (inverse of Ry(B)*Rz(C))
    px = X_neu * cosB * cosC + Y_neu * sinC  - Z_neu * sinB * cosC
    py = -X_neu * cosB * sinC + Y_neu * cosC + Z_neu * sinB * sinC
    pz = X_neu * sinB                        + Z_neu * cosB

    return px, py, pz

def process(input_path, output_path=None):
    out_lines = []
    with open(input_path, "r") as f:
        for raw in f:
            line = raw.rstrip()
            vals = parse_gcode_line(line)
            if vals is None:
                out_lines.append(line)
                continue

            mx, my, mz = vals["X"], vals["Y"], vals["Z"]
            rad_b = math.radians(vals["U"])
            rad_c = math.radians(vals["V"])

            px, py, pz = to_print_position(mx, my, mz, rad_b, rad_c)

            new_line = f"G1 X{px:.4f} Y{py:.4f} Z{pz:.4f} U{vals['U']:.2f} V{vals['V']:.2f}"
            if "E" in vals:
                new_line += f" E{vals['E']:.2f}"
            if "F" in vals:
                new_line += f" F{int(vals['F'])}"
            out_lines.append(new_line)

    result = "\n".join(out_lines)

    if output_path:
        with open(output_path, "w") as f:
            f.write(result)
        print(f"Written to {output_path}")
    else:
        print(result)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tableturn_to_tcp.py <input.txt> [output.txt]")
        sys.exit(1)
    input_file  = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) >= 3 else None
    process(input_file, output_file)
