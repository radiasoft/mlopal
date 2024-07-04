import numpy as np
import sys


def main(data_file):
    data = _read_data(data_file)
    ws, wx = _get_wake(data)
    _write_wake(ws, wx)


def _get_wake(data):
    # TODO(e-carlin): implement by calling model
    a = data["lambda_distribution"]
    return np.random.rand(len(a), len(a[0])), np.random.rand(len(a), len(a[0]))


def _read_data(data_file):
    data = np.fromfile(data_file, dtype=np.float64)
    rows = int(data[0])
    cols = int(data[1])
    num_elements = rows * cols
    return {
        "lambda_distribution": data[2 : num_elements + 2].reshape(rows, cols),
        "bend_radius": data[num_elements + 2],
        "bend_angle": data[num_elements + 3],
        "x_beam_span": data[num_elements + 4],
        "z_beam_span": data[num_elements + 5],
    }


def _write_wake(wake_s, wake_x):
    with open("wake.bin", "wb") as f:
        wake_s.tofile(f)
        wake_x.tofile(f)


if __name__ == "__main__":
    main(sys.argv[1])
