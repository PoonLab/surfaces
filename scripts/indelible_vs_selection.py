import argparse
import csv
import json
import math
import os
from glob import glob


description = """
Compare sitewise HyPhy estimates against INDELible simulated rates and write
one RMSE summary row per simulation.
"""


OMEGAS = [
    0.0510204, 0.1530612, 0.2551020, 0.3571428, 0.4591836, 0.5612244,
    0.6632652, 0.7653060, 0.8673468, 0.9693876, 1.0714284, 1.1734692,
    1.2755100, 1.3775508, 1.4795916, 1.5816324, 1.6836732, 1.7857140,
    1.8877548, 1.9897956, 2.0918364, 2.1938772, 2.2959180, 2.3979588,
    2.4999996, 2.6020404, 2.7040812, 2.8061220, 2.9081628, 3.0102036,
    3.1122444, 3.2142852, 3.3163260, 3.4183668, 3.5204076, 3.6224484,
    3.7244892, 3.8265300, 3.9285708, 4.0306116, 4.1326524, 4.2346932,
    4.3367340, 4.4387748, 4.5408156, 4.6428564, 4.7448972, 4.8469380,
    4.9489788, 5.0510196
]


def parse_args():
    """Command line interface."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "results_dir", type=str,
        help="Path to directory containing method JSON files"
    )
    parser.add_argument(
        "rates_dir", type=str,
        help="Path to directory containing INDELible *_RATES.txt files"
    )
    parser.add_argument(
        "--method", choices=["slac", "fel", "fubar"], required=True,
        help="HyPhy method used to generate the JSON files"
    )
    parser.add_argument(
        "-o", "--outfile", type=argparse.FileType("w"), default=os.sys.stdout,
        help="File to write CSV output (default: stdout)"
    )
    return parser.parse_args()


def load_true_omegas(path):
    """
    Read INDELible *_RATES.txt and map class labels to omega values.
    """
    omegas = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            omega_class = int(row["Class"])
            omegas.append(OMEGAS[omega_class])

    return omegas


def compute_ratio(numerator, denominator):
    """
    Convert method-specific sitewise rates to a comparable dN/dS-style value.
    """
    if numerator is None or denominator is None:
        return None
    if denominator > 0:
        return numerator / denominator
    if denominator == 0 and numerator == 0:
        return 0.0
    return None


def load_slac_sitewise(path):
    """
    Extract sitewise dN/dS values from a HyPhy SLAC JSON file.
    """
    with open(path) as handle:
        data = json.load(handle)

    rows = data["MLE"]["content"]["0"]["by-site"]["AVERAGED"]
    out = []
    for row in rows:
        ds = row[5]
        dn = row[6]
        out.append(compute_ratio(dn, ds))

    return out


def load_fel_sitewise(path):
    """
    Extract sitewise beta/alpha values from a HyPhy FEL JSON file.
    """
    with open(path) as handle:
        data = json.load(handle)

    headers = [header[0] for header in data["MLE"]["headers"]]
    rows = data["MLE"]["content"]["0"]
    alpha_idx = headers.index("alpha")
    beta_idx = headers.index("beta")

    out = []
    for row in rows:
        alpha = row[alpha_idx]
        beta = row[beta_idx]
        out.append(compute_ratio(beta, alpha))

    return out


def load_fubar_sitewise(path):
    """
    Extract sitewise beta/alpha values from a HyPhy FUBAR JSON file.
    """
    with open(path) as handle:
        data = json.load(handle)

    headers = [header[0] for header in data["MLE"]["headers"]]
    rows = data["MLE"]["content"]["0"]
    alpha_idx = headers.index("alpha")
    beta_idx = headers.index("beta")

    out = []
    for row in rows:
        alpha = row[alpha_idx]
        beta = row[beta_idx]
        out.append(compute_ratio(beta, alpha))

    return out


def load_sitewise_estimates(path, method):
    """Dispatch to the parser for the chosen method."""
    if method == "slac":
        return load_slac_sitewise(path)
    if method == "fel":
        return load_fel_sitewise(path)
    return load_fubar_sitewise(path)


def compute_rmse(truth, estimate):
    """Compute true RMSE over non-missing estimated values."""
    sqerr = []
    for true_value, est_value in zip(truth, estimate):
        if est_value is None:
            continue
        sqerr.append((true_value - est_value) ** 2)

    if not sqerr:
        return None, 0

    return math.sqrt(sum(sqerr) / len(sqerr)), len(sqerr)


if __name__ == "__main__":
    args = parse_args()

    suffix = f".{args.method.upper()}.json"
    result_files = sorted(glob(os.path.join(args.results_dir, "*" + suffix)))

    writer = csv.writer(args.outfile)
    writer.writerow([
        "file", "tree.length", "rep", "rmse",
        "n_sites_total", "n_sites_used", "n_sites_dropped"
    ])

    for result_path in result_files:
        filename = os.path.basename(result_path)
        prefix = filename.replace(suffix, "")
        rates_path = os.path.join(args.rates_dir, prefix + "_RATES.txt")

        if not os.path.exists(rates_path):
            os.sys.stderr.write(f"...skipping missing rates file for {prefix}\n")
            continue

        inferred = load_sitewise_estimates(result_path, args.method)
        truth = load_true_omegas(rates_path)

        if len(inferred) != len(truth):
            os.sys.stderr.write(
                f"...skipping length mismatch for {prefix}: "
                f"{len(inferred)} inferred vs {len(truth)} truth\n"
            )
            continue

        rmse, n_used = compute_rmse(truth, inferred)
        if rmse is None:
            os.sys.stderr.write(f"...skipping empty comparison for {prefix}\n")
            continue

        tokens = prefix.split("_")
        tree_length = tokens[1]
        replicate = tokens[2]
        n_total = len(truth)
        n_dropped = n_total - n_used

        writer.writerow([
            prefix, tree_length, replicate, rmse,
            n_total, n_used, n_dropped
        ])
        args.outfile.flush()
