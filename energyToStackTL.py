import numpy as np
import glob
import numpy as np
import warnings
#import pydicom
#import pydicom._storage_sopclass_uids
import os
import datetime as dt
from scipy.signal import convolve2d
import matplotlib.pyplot as plt
from tqdm import tqdm
from PIL import Image
from scipy.ndimage import gaussian_gradient_magnitude
import argparse
import imageio
import SimpleITK as sitk

parser = argparse.ArgumentParser()


parser.add_argument("--base_folder", type=str, default="/home/fismed/Victor/PCI_GAMOS/output/TL/", help="Path to folder that contains Reference/ and Object/ subfolders.")

parser.add_argument("--subfolders", type=str, default="Reference,Object",help="Comma separated subfolder names to process inside base_folder.")

opt = parser.parse_args()

def images_to_stack(input_folder, output_path):
    files_out = []
    stack = []
    files_list = os.listdir(input_folder)
    for file in files_list:

        if "." not in file:
            continue

        ext = file.split(".")[-1].lower()

        if ext == "out":
            files_out.append(file)

        elif ext == "mhd":
            files_out.append(file)

    files_sorted = sorted(files_out, key = lambda file: int(file.split(".")[0].split('_')[1]))
    print(sorted(files_out, key = lambda file: int(file.split(".")[0].split('_')[1])))

    for file in files_sorted:
        try:
            energy_data = np.loadtxt(os.path.join(input_folder, file), delimiter=" ", dtype=np.float32)

        except:
            energy_data_image= sitk.ReadImage(os.path.join(input_folder, file))
            energy_data = sitk.GetArrayFromImage(energy_data_image)

        kernel_size = 5
        sigma = 0.4

        kernel_radius = kernel_size // 2

        # Create a grid of (x,y) coordinates
        x, y = np.meshgrid(np.arange(-kernel_radius, kernel_radius+1), np.arange(-kernel_radius, kernel_radius+1))

        # Compute the Gaussian kernel
        kernel = np.exp(-(x**2 + y**2) / (2 * sigma**2))

        # Normalize the kernel so that its sum is 1
        kernel /= kernel.sum()

        #energy_data = convolve2d(energy_data,kernel, mode = 'same')
        energy_data = energy_data.astype(np.float32)

        stack.append(energy_data)

    stack_array = np.asarray(stack)
    imageio.volwrite(output_path+'Stack.tif', stack_array)
    print(f"Stack saved in: {output_path}")


base = opt.base_folder
if not base.endswith("/"):
    base += "/"

subfolders = [s.strip() for s in opt.subfolders.split(",") if s.strip()]

for sf in subfolders:
    input_folder = os.path.join(base, sf) + "/"
    output_path  = input_folder

    if not os.path.isdir(input_folder):
        print(f"[WARN] Folder not found, skipping: {input_folder}")
        continue

    images_to_stack(input_folder, output_path)

