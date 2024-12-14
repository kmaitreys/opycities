import sys

import numpy as np

from opycities.io import list_materials, usage

material_keys = [
    "ice-w",
    "astrosil",
    "c-gra",
    "c_gra",
    "gra",
    "c-nano",
    "c_nano",
    "c-org",
    "c_org",
    "org",
    "c-p",
    "c_p",
    "c-z",
    "c_z",
    "c",
    "ch3oh-a",
    "ch3oh_a",
    "ch3oh-c",
    "ch3oh_c",
    "ch4-a",
    "ch4_a",
    "ch4-c",
    "ch4_c",
    "co-a",
    "co_a",
    "co",
    "co2-a",
    "co2_a",
    "co2-c",
    "co2_c",
    "co2-w",
    "co2_w",
    "co2",
    "cor-c",
    "cor_c",
    "cor",
    "fe-c",
    "fe_c",
    "iron",
    "fes",
    "tro",
    "h2o-a",
    "h2o_a",
    "h2o-w",
    "h2o_w",
    "h2o",
    "nh3-m",
    "nh3_m",
    "nh3",
    "ol-c-mg00",
    "ol_c_mg00",
    "fay",
    "ol-c-mg100",
    "ol_c_mg100",
    "for",
    "ol-c-mg95",
    "ol_c_mg95",
    "ol-mg40",
    "ol_mg40",
    "ol-mg50",
    "ol_mg50",
    "ol",
    "pyr-c-mg96",
    "pyr_c_mg96",
    "ens",
    "pyr-mg100",
    "pyr_mg100",
    "pyr-mg40",
    "pyr_mg40",
    "pyr-mg50",
    "pyr_mg50",
    "pyr-mg60",
    "pyr_mg60",
    "pyr-mg70",
    "pyr_mg70",
    "pyr",
    "pyr-mg80",
    "pyr_mg80",
    "pyr-mg95",
    "pyr_mg95",
    "sic",
    "sio2",
    "qua",
]


def init_mueller_matrix_dtype(nang):
    return np.dtype(
        [
            ("f11", "f8", nang),
            ("f12", "f8", nang),
            ("f22", "f8", nang),
            ("f33", "f8", nang),
            ("f44", "f8", nang),
            ("f34", "f8", nang),
        ]
    )


def init_particle(nlam, nang):
    fields = np.dtype(
        [
            ("rv", "f8"),
            ("rvmin", "f8"),
            ("rvmax", "f8"),
            ("rho", "f8"),
            ("Kabs", "f8", nlam),
            ("Ksca", "f8", nlam),
            ("Kext", "f8", nlam),
            ("g", "f8", nlam),
            ("F", init_mueller_matrix_dtype(nang), nlam),
            ("testscat", "bool", nlam),
            ("scat_ok", "bool"),
            ("scat_ok_lmin", "f8"),
        ]
    )
    return np.zeros((1,), dtype=fields)


def add_material(material, mfrac=None):
    print(f"Adding material '{material}' with mass fraction: {mfrac}")


def set_porosity(porosity, pmantle=None):
    print(f"Setting porosity: {porosity} (Mantle Porosity: {pmantle})")


def output_fits(output_directory):
    print(f"Outputting to FITS in directory: {output_directory}")


def main():
    mat_loc = [None for i in range(21)]
    mat_lnk = [None for i in range(21)]
    mat_cmd = [None for i in range(21)]
    mat_rfn = [None for i in range(21)]
    mat_rfk = [None for i in range(21)]
    mat_rho = [None for i in range(21)]
    mat_mfr = [None for i in range(21)]

    amin = 0.05
    amax = 3000
    apow = 3.5
    amean = 0
    asig = 0
    na = 0
    sdkind = "apow"
    lmin = 0.05
    lmax = 10000.0
    nlam = 300
    nang = 180
    chopangle = 0.0
    pcore = 0
    pmantle = 0
    nm = 0
    nmant = 0
    method = "DHS"
    fmax = 0.8
    write_fits = False
    write_scatter = False
    for_radmc = False

    mat_rho = np.zeros(21, dtype=np.float64)

    # Parse command-line arguments
    for arg in sys.argv[1:]:
        if arg == "-q":
            quiet = True
        if arg == "-v":
            verbose = True
        if arg == "-debug":
            debug = True

        if arg.startswith("--"):
            # trim off the leading '--'
            arg = arg[2:]
        match arg:
            case "help" | "-h" | "-help":
                # See there are no more arguments
                if len(sys.argv) > 2:
                    print("Error: --help must be the only argument")
                    sys.exit(1)
                else:
                    usage()
                    sys.exit(0)

            case "-c" | "-m":
                nm += 1
                if nm > 20:
                    print("Too many materials")
                    sys.exit(1)

                if len(sys.argv) <= 2:
                    list_materials()
                    sys.exit(0)

                if arg == "-m":
                    mat_loc[nm] = "mantle"
                    nmant += 1
                else:
                    mat_loc[nm] = "core"

                # See if next arg is 3 numbers separated by colons
                if ":" in sys.argv[sys.argv.index(arg) + 1]:
                    try:
                        n, k, rho = sys.argv[sys.argv.index(arg) + 1].split(":")
                    except ValueError:
                        print(
                            "Error: -c and -m arguments must be followed by 3 numbers separated by colons"
                        )
                        sys.exit(1)
                    print(f"n: {n}, k: {k}, rho: {rho}")

                    mat_rfn[nm] = float(n)
                    mat_rfk[nm] = float(k)
                    mat_rho[nm] = float(rho)
                    if mat_rho[nm] < 0.0:
                        print("Error: Density must be positive")
                        sys.exit(1)

                    mat_cmd[nm] = True
                else:
                    mat_cmd[nm] = False
                    # See if next arg is a material key
                    if sys.argv[sys.argv.index(arg) + 1] not in material_keys:
                        print("Error: Material key not recognized")
                        sys.exit(1)

                # Next value will be mass fraction
                # Ensure it is a float
                try:
                    mfrac = float(sys.argv[sys.argv.index(arg) + 2])
                    mat_mfr[nm] = mfrac
                except ValueError:
                    print("No mass fraction given. Defaulting to 1.0")
                    mat_mfr[nm] = 1.0
                    # goto next arg
                    continue

                if mat_mfr[nm] == 0.0:
                    print("Warning: Ignoring material with mass fraction of 0")
                    nm -= 1

                # There might be a next arg (density)
                if sys.argv[sys.argv.index(arg) + 3]:
                    try:
                        mat_rho[nm] = float(sys.argv[sys.argv.index(arg) + 3])
                    except ValueError:
                        print("Error: Density must be a float")
                        sys.exit(1)

                    if mat_cmd[nm]:
                        print("Error: Cannot specify density twice")
                        sys.exit(1)

                    if mat_rho[nm] < 0.0:
                        print("Error: Density must be positive")
                        sys.exit(1)

            case "-diana" | "-dsharp" | "-dsharp-no-ice":
                ...


if __name__ == "__main__":
    main()
