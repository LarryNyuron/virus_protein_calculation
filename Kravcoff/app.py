# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, html, dcc
import plotly.express as px
import pandas as pd

app = Dash(__name__)

# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options
df = pd.DataFrame({
    "X": [244.125, 243.246, 242.478, 241.898, 244.092, 243.283, 244.159, 243.66, 245.474, 242.481, 241.8, 240.757, 239.933, 242.828, 244.017, 243.824, 245.143, 240.795, 239.837, 
240.314, 241.317, 238.486, 238.559, 237.205, 236.887, 236.489, 236.346, 236.23, 239.587, 239.865, 240.239, 241.348, 240.964, 241.341, 242.19, 239.3, 239.495, 240.655, 241.575, 238.164, 237.859, 238.241, 240.605, 241.623, 243.058, 243.911, 241.463, 241.066, 240.43, 251.257, 251.944, 253.27, 253.525, 252.188, 252.818, 254.129, 255.419, 255.219, 254.434, 256.055, 255.437, 253.999, 255.924, 255.802, 256.254, 257.362, 256.688, 256.644, 257.474, 257.186, 255.395, 255.724, 255.652, 255.082, 254.713, 253.435, 254.588, 256.241, 256.221, 255.811, 255.603, 
257.615, 258.694, 257.624, 258.804, 255.692, 255.338, 256.223, 256.21, 253.868, 252.951, 251.511, 250.623, 250.735, 257.016, 257.917, 257.515, 257.496, 259.345, 260.041, 260.001, 259.889, 260.103, 257.168, 256.763, 258.028, 259.1, 256.021, 256.833, 257.069, 257.921, 259.106, 258.84, 257.767, 259.981, 259.454, 259.425, 258.964, 258.91, 258.446, 258.42, 259.853, 259.808, 261.228, 262.117, 259.279, 259.509, 258.843, 258.981, 258.113, 261.453, 262.799, 262.769, 261.871, 263.475, 262.839, 263.745, 263.857, 264.275, 265.352, 264.898, 264.624, 265.645, 
266.854, 265.241, 263.405, 263.663, 264.259, 265.092, 262.35, 261.975, 262.493, 263.856, 264.34, 264.244, 263.248, 263.539, 263.991, 265.327, 262.956, 265.285, 265.325, 265.301, 266.261, 266.631, 266.575, 265.593, 267.561, 265.603, 267.579, 266.595, 264.247, 264.231, 263.603, 263.151, 263.593, 263.594, 263.029, 264.05, 263.715, 265.298, 266.392, 266.864, 267.714, 267.596, 267.878, 267.321, 266.336, 266.74, 265.821, 265.584, 268.173, 268.422, 269.703, 269.528, 269.097, 265.287, 264.404, 265.057, 266.267, 264.207, 263.335, 263.564, 262.422, 264.271, 264.811, 263.751, 262.685, 265.713, 264.068, 263.19, 263.775, 264.984, 263.117, 262.511, 262.126, 262.977, 260.973, 262.915, 263.343, 262.474, 261.249, 263.221, 263.102, 262.341, 262.889, 263.612, 262.404, 261.729, 262.536, 263.018, 262.004, 260.808, 263.5, 262.573, 261.376, 263.341, 262.499, 261.649, 261.723, 262.639, 262.08, 261.781, 263.567, 260.739, 260.674, 261.066, 261.284, 259.269, 258.264, 258.862, 258.344, 261.133, 261.529, 262.987, 263.415, 261.362, 263.736, 265.152, 265.604, 266.793, 266.025, 266.124, 266.759, 265.498, 264.683, 265.061, 265.65, 265.174, 263.863, 262.693, 263.47, 261.436, 266.699, 267.413, 266.645, 266.97, 267.836, 266.668, 266.049, 266.371, 265.634, 264.809, 264.461, 264.926, 263.533, 262.834, 263.645, 263.262, 264.491, 264.504, 262.632, 262.5, 261.276, 265.524, 266.76, 267.222, 267.306, 267.856, 269.089, 267.507, 267.974, 266.966, 267.344, 268.316, 267.107, 269.116, 265.684, 264.673, 264.545, 264.552, 263.314, 263.12, 261.681, 263.472, 264.419, 264.264, 265.558, 265.574, 263.857, 264.983, 262.689, 266.638, 267.957, 268.05, 268.878, 268.917, 268.995, 270.312, 267.224, 267.283, 266.167, 265.932, 267.332, 268.515, 269.767, 268.381, 270.874, 269.476, 270.732, 265.489, 264.399, 264.555, 265.359, 263.061, 262.855, 262.709, 262.933, 262.662, 262.886, 262.76, 262.804, 263.804, 263.889, 262.615, 262.561, 264.107, 265.283, 265.871, 264.841, 265.092, 266.343, 264.094, 261.588, 260.306, 259.846, 259.78, 259.26, 259.661, 260.036, 
259.716, 260.304, 260.117, 259.515, 259.068, 258.024, 258.234, 260.253, 256.901, 255.824, 255.352, 255.448, 254.641, 253.591, 254.836, 254.362, 252.93, 252.031, 254.442, 255.83, 255.688, 256.506, 252.722, 251.389, 250.74, 249.599, 251.458, 251.756, 250.704, 249.536, 251.042, 251.483, 251.001, 252.156, 253.312, 250.606, 251.727, 251.862, 252.931, 252.404, 251.343, 253.466, 254.511, 255.85, 254.575, 253.14, 252.689, 253.775, 254.664, 251.483, 251.838, 252.133, 251.919, 252.39, 252.266, 251.735, 252.425, 251.894, 252.238, 253.687, 254.673, 254.143, 
252.932, 255.105, 253.93, 256.263, 255.064, 254.717, 255.752, 256.951, 254.526, 254.702, 255.504, 255.264, 256.109, 255.908, 254.824, 255.796, 254.394, 256.819, 254.154, 256.965, 256.939, 257.352, 258.545, 257.887, 257.63, 257.585, 257.432, 257.373, 257.276, 256.372, 256.695, 257.543, 257.35, 255.326, 254.552, 254.943, 258.486, 259.348, 259.405, 258.689, 260.782, 261.385, 260.783, 260.289, 260.505, 260.561, 259.964, 261.788, 262.026, 260.727, 263.109, 261.283, 261.446, 260.165, 259.864, 262.575, 263.868, 264.45, 265.406, 263.881, 259.421, 258.18, 257.427, 257.089, 257.316, 257.174, 256.45, 257.134, 256.477, 256.283, 256.207, 257.344, 258.453, 259.165, 259.517, 260.048, 260.456, 261.337, 261.155, 259.203, 259.537, 258.579, 258.866, 260.93, 262.001, 262.514, 262.499, 263.511, 263.494, 264.004, 257.426, 256.473, 257.248, 258.236, 255.538, 255.476, 256.926, 256.829, 257.537, 256.633, 255.945, 258.615, 259.525, 
259.374, 256.647, 255.827, 256.508, 257.715, 255.636, 254.865, 254.923, 255.733, 256.253, 255.196, 254.015, 256.68, 255.538, 257.122, 255.632, 254.725, 255.151, 256.345, 254.157, 254.38, 253.478, 252.493, 254.053, 255.188, 252.735, 253.815, 253.02, 253.477, 254.625, 253.168, 252.161, 252.588, 252.923, 252.67, 251.64, 252.083, 252.553, 251.875, 253.768, 252.589, 253.76, 254.866, 254.8, 255.905, 255.862, 253.611, 253.532, 254.249, 255.258, 254.156, 254.567, 253.177, 253.731, 254.361, 255.756, 255.884, 253.458, 252.122, 252.453, 256.792, 258.17, 258.344, 259.266, 259.1, 257.462, 257.533, 256.695, 256.088, 257.03, 255.601, 254.858, 253.753, 255.457, 256.671, 255.882, 256.7, 257.59, 254.832, 254.142, 256.369, 257.045, 256.955, 257.641, 256.472, 255.095, 255.329, 256.143, 256.038, 257.285, 258.012, 254.735, 254.363, 254.919, 257.525, 258.654, 258.128, 257.079, 259.546, 259.339, 259.201, 258.877, 258.471, 258.008, 257.074, 259.719, 260.327, 260.216, 258.654, 258.28, 256.853, 256.359, 259.223, 256.189, 254.815, 253.934, 253.28, 254.818, 255.394, 255.751, 256.105, 255.666, 253.902, 253.096, 251.602, 250.815, 253.548, 254.858, 252.478, 255.341, 251.208, 249.795, 249.208, 248.017, 249.553, 250.651, 249.389, 250.044, 249.559, 249.238, 248.56, 250.622, 250.948, 251.967, 252.435, 253.495, 249.724, 249.517, 248.35, 247.962, 250.793, 251.854, 250.602, 247.796, 246.664, 247.03, 246.249, 246.095, 244.884, 243.728, 244.899, 242.614, 243.788, 242.649, 241.547, 248.207, 248.656, 249.041, 249.494, 248.859, 249.18, 248.816, 247.867, 249.566, 249.336, 250.379, 251.577, 249.47, 249.177, 247.722, 247.204, 247.051, 249.932, 250.856, 250.533, 249.383, 250.836, 251.938, 250.987, 251.986, 251.553, 251.357, 252.078, 252.781, 251.839, 251.058, 249.722, 251.663, 249.01, 250.955, 249.625, 251.902, 252.553, 253.275, 253.356, 251.527, 252.294, 
253.8, 254.525, 255.027, 255.699, 255.748, 255.294, 256.652, 254.608, 254.707, 255.128, 254.482, 255.152, 253.17, 252.456, 252.535, 251.59, 253.68, 253.961, 252.891, 251.699, 254.148, 
253.329, 252.404, 251.299, 250.249, 251.838, 252.494, 250.324, 253.979, 251.585, 250.646, 250.584, 249.683, 250.945, 252.437, 253.227, 252.97, 254.195, 251.545, 251.671, 250.579, 249.542, 251.578, 252.423, 250.153, 250.834, 249.935, 250.791, 251.975, 248.924, 247.525, 246.72, 246.806, 250.191, 250.905, 251.453, 251.184, 249.955, 249.315, 252.253, 252.808, 251.754, 250.851, 253.941, 254.345, 253.012, 251.871, 250.95, 251.243, 252.4, 251.142, 250.697, 251.464, 249.203, 250.194, 250.357, 249.808, 248.881, 249.605, 250.303, 249.54, 251.734, 250.381, 
249.892, 249.978, 251.014, 250.678, 250.335, 250.327, 248.874, 248.812, 249.109, 248.548, 247.426, 247.051, 248.028, 247.76, 248.809, 250.013, 250.359, 249.327, 249.064, 251.72, 251.97, 248.726, 247.724, 248.354, 248.568, 247.103, 248.26, 248.902, 248.641, 249.262, 248.418, 248.945, 249.617, 250.91, 251.246, 252.04, 247.106, 246.228, 246.246, 245.231, 244.81, 244.181, 244.979, 245.188, 245.397, 247.42, 247.588, 249.036, 249.442, 247.174, 246.805, 245.132, 245.395, 249.803, 251.2, 251.828, 251.612, 251.972, 251.828, 252.389, 250.888, 252.62, 253.276, 253.938, 255.015, 254.318, 254.329, 254.742, 254.159, 255.75, 253.301, 253.812, 254.885, 255.32, 252.565, 251.83, 251.964, 255.29, 256.359, 257.65, 258.05, 256.36, 255.173, 255.458, 255.227, 254.029, 252.944, 253.91, 258.3, 259.518, 260.749, 261.839, 259.658, 258.432, 259.835, 260.591, 261.718, 261.225, 260.262, 262.474, 263.625, 264.579, 265.849, 266.882, 261.874, 261.442, 262.488, 263.651, 260.186, 259.373, 258.92, 259.19, 262.052, 262.922, 263.176, 264.116, 262.29, 262.172, 262.326, 262.455, 262.517, 262.238, 261.237, 259.972, 261.262, 258.729, 262.875, 262.954, 261.654, 261.361, 264.127, 264.312, 264.562, 265.693, 263.499, 260.865, 259.619, 259.83, 260.847, 258.563, 258.098, 257.353, 258.39, 256.908, 257.952, 257.212, 
256.764, 258.866, 258.996, 257.991, 258.363, 258.832, 259.226, 258.945, 260.687, 256.719, 255.644, 255.641, 254.589, 254.283, 254.309, 254.861, 253.763, 256.818, 256.977, 256.745, 256.867, 258.384, 258.676, 256.406, 256.17, 257.396, 258.478, 255.639, 254.944, 255.93, 257.2, 258.223, 257.656, 256.598, 258.715, 259.76, 260.405, 259.451, 260.191, 258.371, 257.925, 258.481, 259.632, 258.277, 257.918, 256.404, 258.438, 257.65, 258.06, 257.909, 256.808, 257.19, 257.269, 256.433, 258.727, 259.011, 258.988, 259.326, 260.046, 259.996, 259.667, 259.937, 260.658, 258.807, 259.076, 259.088, 258.065, 258.029, 258.263, 260.253, 260.408, 260.794, 261.578, 261.506, 260.967, 262.012, 260.789, 260.234, 260.545, 260.654, 259.681, 259.443, 258.178, 259.386, 261.852, 262.113, 261.809, 261.687, 263.561, 261.678, 261.418, 262.487, 263.605, 261.542, 260.231, 259.503, 260.11, 258.196, 262.166, 263.205, 264.427, 265.571, 262.494, 261.469, 260.938, 264.163, 265.204, 265.402, 264.787, 264.806, 264.318, 266.004, 266.268, 266.57, 266.429, 266.7, 267.971, 265.984, 265.793, 267.094, 267.987, 264.703, 264.096, 265.285, 
267.219, 268.434, 268.654, 269.7, 268.116, 267.194, 266.282, 267.645, 267.68, 266.5, 265.835, 267.612, 266.241, 265.143, 263.789, 263.65, 265.339, 266.538, 262.791, 261.452, 260.659, 259.614, 260.652, 260.987, 260.94, 261.142, 260.412, 261.237, 262.289, 259.163, 259.538, 260.747, 261.423, 260.431, 259.337, 261.851, 262.559, 262.744, 262.915, 260.807, 259.936, 260.354, 261.454, 260.043, 259.723, 259.093, 259.923, 259.481, 259.798, 258.822, 257.715, 259.794, 258.515, 260.864, 259.248, 258.412, 258.208, 259.173, 259.069, 258.423, 258.912, 256.943, 256.584, 255.757, 254.993, 255.776, 256.498, 255.904, 255.155, 255.375, 256.088, 254.759, 254.876, 255.249, 254.706, 253.551, 252.943, 253.775, 256.183, 256.598, 256.126, 256.082, 258.114, 258.695, 260.2, 258.136, 255.755, 255.322, 256.414, 256.524, 253.973, 253.465, 257.233, 258.332, 258.284, 257.452, 259.632, 259.628, 261.036, 260.463, 259.175, 259.265, 260.228, 260.464, 257.922, 257.559, 256.858, 257.895, 256.778, 257.4, 260.714, 261.69, 263.06, 263.477, 261.359, 260.647, 263.776, 265.099, 266.087, 266.116, 265.508, 264.921, 263.521, 266.894, 267.923, 269.232, 269.629, 267.976, 269.221, 269.512, 269.01, 269.903, 271.153, 272.315, 273.171, 271.552, 270.395, 272.787, 269.877, 272.342, 273.402, 274.211, 273.996, 272.815, 273.873, 272.063, 275.158, 275.969, 275.222, 274.828, 277.306, 277.993, 278.276, 278.247, 275.012, 274.326, 275.324, 274.957, 273.366, 273.976, 272.067, 276.588, 277.635, 277.111, 276.805, 278.893, 279.901, 277.004, 276.507, 275.217, 274.847, 277.519, 277.846, 278.79, 274.6],
    "Y": [614.17, 615.326, 615.763, 616.854, 616.487, 617.646, 618.716, 619.722, 618.507, 614.892, 615.148, 614.066, 613.75, 615.197, 616.071, 617.289, 615.533, 613.495, 612.457, 611.563, 610.859, 613.089, 614.205, 614.447, 613.468, 612.213, 611.749, 611.417, 611.607, 610.818, 609.358, 609.048, 611.475, 610.566, 611.842, 608.465, 607.029, 606.495, 607.244, 606.264, 606.349, 
604.814, 605.203, 604.54, 604.833, 603.941, 604.901, 606.272, 604.012, 586.741, 585.678, 586.165, 585.99, 584.479, 583.407, 586.789, 587.283, 588.33, 589.262, 587.882, 587.099, 587.043, 588.177, 589.135, 590.491, 590.625, 588.751, 587.336, 586.404, 587.355, 591.492, 592.845, 593.708, 593.291, 593.394, 593.542, 592.453, 594.897, 595.85, 597.195, 597.338, 596.038, 595.423, 595.469, 596.072, 598.175, 599.514, 600.481, 600.519, 599.788, 599.161, 599.55, 599.092, 597.624, 601.24, 602.219, 603.611, 603.91, 601.942, 600.905, 601.27, 602.448, 600.269, 604.479, 605.859, 606.678, 606.251, 606.249, 605.556, 604.183, 607.844, 608.668, 610.107, 610.436, 608.024, 608.205, 609.469, 607.126, 609.658, 607.302, 608.576, 610.948, 612.347, 612.757, 612.481, 613.254, 614.726, 615.677, 616.895, 615.13, 613.404, 613.849, 615.242, 615.564, 612.859, 612.844, 616.066, 617.428, 617.331, 616.841, 618.207, 618.221, 619.023, 618.794, 619.878, 617.799, 617.74, 619.036, 619.021, 617.403, 618.504, 616.062, 620.155, 621.443, 622.555, 622.675, 621.852, 623.124, 622.849, 623.558, 623.374, 624.509, 625.731, 625.956, 624.522, 623.83, 624.155, 622.914, 623.589, 622.342, 622.682, 626.539, 627.686, 628.984, 629.111, 627.275, 629.949, 631.256, 632.347, 633.535, 631.932, 632.867, 633.057, 633.905, 632.453, 631.071, 632.67, 632.273, 632.414, 631.625, 630.446, 631.928, 630.542, 629.944, 629.483, 630.577, 632.275, 631.589, 630.26, 630.209, 632.433, 633.639, 634.34, 633.888, 629.184, 627.869, 626.817, 626.818, 627.385, 625.924, 624.824, 623.553, 623.459, 624.662, 625.828, 625.454, 624.885, 625.727, 622.579, 621.308, 620.198, 620.276, 621.348, 619.164, 618.06, 616.737, 616.053, 618.102, 616.996, 616.371, 615.126, 613.992, 614.216, 615.398, 616.325, 615.559, 616.929, 612.77, 611.586, 610.907, 610.125, 610.585, 609.26, 610.723, 611.206, 610.682, 609.227, 608.708, 610.872, 609.98, 612.339, 608.534, 608.562, 607.172, 607.039, 605.949, 606.603, 608.147, 608.146, 609.098, 609.407, 608.402, 607.207, 607.252, 606.114, 609.564, 610.469, 609.732, 608.665, 611.31, 610.397, 612.274, 611.143, 610.321, 609.792, 608.756, 607.573, 610.954, 611.867, 612.466, 611.989, 609.215, 608.346, 607.08, 605.994, 609.078, 609.529, 607.223, 606.082, 605.236, 604.028, 606.541, 605.357, 607.168, 605.88, 605.188, 604.281, 603.062, 606.203, 605.554, 604.877, 604.123, 603.098, 602.029, 605.054, 605.521, 606.242, 603.414, 602.485, 601.228, 600.125, 603.16, 603.964, 604.451, 603.085, 601.394, 600.256, 599.459, 598.357, 600.72, 601.352, 601.705, 600.028, 599.413, 597.92, 597.223, 600.139, 601.524, 599.512, 597.425, 596.012, 595.129, 594.037, 595.873, 596.521, 595.943, 597.713, 596.54, 598.321, 597.732, 595.588, 594.823, 594.788, 595.527, 595.463, 595.513, 594.339, 596.713, 594.349, 596.739, 595.546, 595.53, 593.923, 593.859, 594.398, 594.667, 592.43, 591.744, 590.705, 589.877, 588.784, 588.389, 588.077, 594.547, 595.062, 596.161, 595.953, 593.958, 592.746, 592.803, 591.442, 591.583, 590.738, 597.322, 598.443, 599.259, 599.654, 599.329, 599.51, 600.305, 601.414, 601.291, 599.414, 600.19, 602.496, 603.612, 603.421, 603.166, 604.903, 605.451, 606.575, 605.947, 603.544, 603.418, 604.785, 604.922, 602.92, 601.446, 600.612, 600.654, 599.923, 605.797, 607.173, 608.037, 607.67, 607.58, 607.523, 609.184, 610.059, 611.449, 611.626, 609.573, 610.445, 610.176, 610.13, 612.442, 613.807, 614.807, 614.532, 614.067, 614.224, 615.376, 613.196, 615.138, 613.8, 611.818, 613.078, 611.091, 611.726, 615.977, 617.021, 618.09, 618.329, 617.672, 618.39, 618.617, 618.714, 619.767, 620.879, 620.647, 619.225, 620.299, 618.1, 622.092, 623.272, 623.946, 624.443, 624.285, 624.045, 624.16, 624.73, 623.948, 624.551, 626.013, 626.324, 623.811, 622.344, 621.541, 621.526, 620.291, 620.255, 626.93, 628.357, 628.615, 627.967, 629.017, 628.044, 626.722, 629.542, 629.834, 631.304, 632.13, 629.367, 630.192, 627.916, 631.609, 632.954, 634.061, 635.127, 632.952, 634.221, 634.655, 633.966, 633.798, 634.747, 634.879, 635.935, 634.286, 633.966, 635.16, 635.024, 636.339, 633.788, 633.773, 635.097, 635.728, 632.613, 635.549, 636.812, 638.039, 638.928, 636.922, 635.489, 634.894, 638.094, 639.248, 639.144, 640.084, 639.516, 640.357, 638.205, 638.01, 637.821, 636.913, 635.729, 637.225, 638.108, 639.133, 637.917, 639.954, 638.732, 639.753, 637.442, 636.594, 635.795, 636.291, 637.607, 638.664, 638.821, 634.563, 633.747, 632.804, 631.978, 632.878, 633.702, 632.11, 632.911, 632.029, 630.668, 630.577, 632.59, 633.795, 631.577, 629.612, 628.253, 627.329, 627.389, 627.747, 627.901, 626.301, 626.494, 625.577, 624.124, 623.804, 623.243, 621.807, 621.127, 621.693, 621.23, 621.528, 621.81, 619.896, 619.164, 617.72, 617.43, 619.795, 619.0, 616.813, 615.404, 614.874, 615.163, 614.668, 614.881, 615.491, 614.392, 615.399, 614.727, 613.695, 614.384, 613.351, 613.697, 614.082, 613.562, 612.222, 611.946, 614.602, 613.951, 615.735, 611.37, 610.071, 610.338, 610.976, 609.449, 610.044, 611.485, 609.864, 610.072, 609.926, 610.499, 
609.105, 609.157, 608.946, 609.979, 609.673, 607.552, 607.308, 606.469, 605.984, 606.299, 611.203, 612.269, 613.444, 613.909, 612.78, 613.887, 613.969, 615.105, 616.486, 617.407, 615.1, 614.572, 613.427, 616.654, 617.978, 618.454, 617.657, 618.074, 616.713, 619.048, 619.767, 620.421, 621.485, 622.073, 621.145, 622.558, 620.719, 621.783, 622.781, 624.11, 624.724, 622.936, 621.583, 621.264, 624.56, 625.831, 625.79, 626.777, 626.169, 624.651, 624.49, 623.699, 622.736, 623.806, 622.415, 621.943, 620.781, 622.85, 624.111, 623.41, 623.681, 622.74, 623.754, 623.029, 623.356, 623.209, 624.954, 625.31, 624.839, 624.988, 626.831, 627.535, 627.241, 624.286, 623.769, 622.291, 621.662, 623.893, 625.302, 625.251, 626.631, 626.538, 621.746, 
620.339, 620.111, 621.022, 619.764, 619.87, 618.315, 618.897, 618.551, 618.559, 619.018, 617.175, 616.797, 617.571, 615.673, 617.24, 615.331, 616.12, 615.79, 618.031, 618.026, 619.448, 619.704, 620.367, 621.772, 622.25, 621.755, 623.221, 623.742, 624.798, 624.545, 622.612, 622.994, 623.328, 623.222, 623.744, 625.978, 627.055, 627.694, 627.701, 628.132, 627.851, 629.508, 628.859, 628.244, 628.874, 630.191, 630.641, 627.954, 626.693, 626.724, 625.469, 625.558, 624.294, 624.338, 630.8, 632.054, 631.945, 632.919, 633.174, 634.789, 630.75, 630.477, 631.765, 632.56, 629.558, 628.121, 629.607, 627.473, 631.966, 633.186, 634.317, 635.235, 634.226, 635.28, 636.457, 637.212, 636.577, 637.647, 637.806, 637.797, 638.957, 637.93, 638.122, 
637.05, 637.186, 639.573, 640.443, 639.58, 640.682, 635.952, 634.845, 633.847, 633.915, 634.098, 633.977, 633.112, 631.913, 633.721, 632.924, 631.853, 630.82, 631.113, 632.357, 631.548, 632.235, 629.62, 628.466, 627.23, 627.324, 628.563, 629.133, 629.242, 628.238, 626.07, 624.818, 624.765, 625.649, 623.637, 623.73, 623.737, 623.623, 623.146, 622.376, 622.611, 622.8, 622.905, 623.622, 623.24, 621.808, 621.402, 624.141, 625.587, 626.479, 625.682, 621.043, 619.657, 619.447, 620.139, 618.728, 618.705, 617.342, 618.212, 618.49, 618.188, 616.688, 616.163, 618.958, 618.432, 620.432, 616.004, 614.566, 614.217, 614.81, 614.034, 614.245, 613.539, 613.887, 613.389, 613.267, 612.825, 611.808, 610.847, 612.156, 611.202, 611.997, 611.04, 609.7, 609.427, 611.765, 612.496, 613.082, 608.868, 607.577, 606.643, 605.823, 606.899, 607.427, 606.634, 607.322, 606.746, 605.894, 606.469, 606.494, 605.927, 607.318, 608.308, 608.054, 609.347, 606.931, 607.55, 607.434, 608.159, 609.029, 609.72, 610.424, 611.889, 606.513, 606.332, 605.184, 605.06, 607.607, 608.066, 609.749, 610.65, 604.362, 603.206, 603.533, 604.13, 602.64, 601.135, 600.505, 600.775, 599.648, 603.129, 603.377, 602.41, 602.477, 603.285, 602.193, 602.516, 601.494, 600.569, 601.354, 601.647, 599.339, 598.426, 597.028, 596.891, 596.959, 597.167, 596.803, 601.717, 602.479, 601.6, 602.102, 603.434, 604.309, 602.657, 600.299, 599.368, 597.986, 597.478, 599.262, 598.269, 598.426, 597.597, 597.93, 597.371, 596.046, 595.243, 595.624, 596.172, 594.91, 594.508, 594.322, 594.121, 593.237, 593.811, 593.418, 591.848, 591.203, 594.742, 595.373, 596.864, 597.337, 595.074, 595.363, 593.629, 594.893, 597.607, 599.058, 599.694, 599.621, 599.561, 601.046, 601.477, 601.427, 601.876, 600.306, 600.918, 602.351, 602.958, 600.905, 599.536, 598.771, 599.01, 597.511, 597.761, 597.015, 595.782, 602.901, 604.265, 605.223, 606.302, 604.297, 603.011, 603.17, 602.703, 604.83, 605.667, 605.887, 606.138, 605.113, 603.61, 602.862, 603.178, 605.804, 606.012, 607.483, 608.355, 605.613, 606.148, 607.781, 609.175, 610.084, 609.732, 609.044, 607.731, 606.861, 611.244, 612.278, 613.647, 614.028, 612.352, 613.429, 613.422, 613.898, 614.097, 614.395, 615.721, 616.91, 616.913, 615.905, 617.236, 617.421, 617.255, 617.927, 619.129, 620.324, 620.599, 619.357, 618.272, 618.658, 618.078, 621.033, 622.216, 623.443, 623.354, 622.102, 620.888, 623.334, 620.635, 624.592, 625.79, 627.037, 627.71, 625.927, 627.073, 627.342, 628.514, 629.739, 629.65, 628.307, 627.494, 629.64, 626.06, 630.881, 632.153, 633.148, 633.797, 632.591, 632.586, 631.65, 633.244, 634.113, 635.585, 636.069, 633.953, 636.282, 637.708, 638.343, 637.837, 638.24, 638.494, 639.655, 640.679, 639.514, 639.45, 640.063, 640.438, 640.305, 641.279, 641.62, 640.261, 640.887, 641.285, 640.17, 640.176, 642.556, 643.536, 643.114, 639.223, 638.061, 638.24, 639.309, 637.551, 637.181, 637.159, 636.715, 636.194, 636.105, 635.93, 636.045, 636.908, 636.505, 635.007, 634.566, 636.869, 638.041, 637.586, 634.236, 632.785, 632.322, 633.137, 632.153, 631.016, 630.454, 630.725, 630.717, 628.948, 628.672, 630.962, 631.245, 629.991, 630.063, 632.047, 631.589, 633.529, 628.838, 627.612, 626.372, 626.429, 627.576, 627.803, 625.244, 623.976, 623.035, 622.818, 623.37, 624.428, 622.179, 623.972, 622.484, 621.572, 620.144, 619.733, 621.8, 623.245, 620.906, 623.574, 619.383, 618.0, 617.03, 617.401, 617.764, 618.119, 618.607, 615.781, 614.748, 613.671, 613.154, 614.15, 612.838, 615.109, 613.36, 612.355, 611.246, 611.47, 612.989, 614.026, 610.047, 608.92, 607.667, 607.681, 606.576, 605.297, 604.26, 604.251, 604.897, 606.056, 603.853, 603.388, 602.365, 600.999, 600.686, 602.343, 603.518, 603.4, 603.528, 600.192, 598.834, 597.997, 597.954, 598.57, 597.307, 597.343, 596.562, 595.108, 594.715, 597.172, 598.659, 599.396, 599.488, 594.314, 592.887, 592.176, 592.63, 592.211, 591.936, 592.834, 590.913, 592.375, 591.212, 591.02, 590.179, 590.81, 591.684, 590.027, 588.817, 590.363, 590.871, 590.786, 589.804, 589.995, 588.686, 589.116, 591.834, 591.912, 591.59, 592.287, 593.316, 593.617, 592.465, 594.927, 590.536, 590.137, 591.095, 590.867, 588.719, 587.767, 588.258, 587.849, 592.173, 593.157, 593.161, 592.33, 594.561, 595.474, 595.068, 594.082, 594.196, 595.044, 596.174, 594.856, 595.419, 594.64, 596.644, 594.499, 595.244, 595.846, 596.249, 594.351, 593.073, 594.185, 595.894, 596.461, 597.772, 597.882, 596.72, 597.366, 598.765, 600.07, 600.322, 599.45, 601.183, 601.121, 601.023, 601.384],
    "Z": [160.997, 161.339, 160.097, 160.055, 161.862, 162.409, 163.031, 163.544, 162.989, 159.09, 157.825, 157.545, 158.407, 156.694, 157.037, 157.248, 157.106, 156.344, 155.985, 154.842, 154.982, 155.624, 154.591, 153.925, 152.879, 153.091, 154.325, 152.06, 153.722, 152.516, 152.773, 153.215, 151.629, 150.586, 152.447, 152.476, 152.664, 151.799, 151.445, 152.34, 150.852, 152.814, 151.473, 150.658, 151.093, 151.092, 149.154, 149.014, 148.508, 165.563, 166.288, 166.879, 168.073, 165.361, 166.047, 166.053, 166.542, 167.628, 167.451, 165.291, 164.182, 164.61, 168.746, 169.837, 169.334, 168.804, 171.017, 171.574, 170.691, 172.999, 169.495, 169.074, 170.312, 171.318, 168.053, 168.683, 166.864, 170.239, 171.348, 170.778, 169.571, 172.023, 171.146, 173.425, 169.771, 171.66, 171.254, 172.012, 173.243, 171.541, 170.525, 170.774, 169.63, 169.369, 171.264, 171.84, 171.401, 170.208, 171.392, 172.228, 173.69, 174.044, 174.552, 
172.362, 172.083, 171.889, 172.326, 173.346, 174.398, 173.801, 171.254, 171.029, 170.6, 170.093, 169.96, 168.57, 167.979, 167.858, 166.703, 166.578, 165.998, 170.78, 170.376, 170.082, 
170.881, 171.479, 171.151, 172.123, 172.007, 173.085, 168.943, 168.598, 167.977, 167.193, 167.644, 166.382, 168.359, 167.858, 166.414, 166.093, 168.649, 170.13, 170.891, 170.681, 171.704, 165.54, 164.122, 163.598, 162.691, 163.382, 162.4, 162.697, 164.188, 163.728, 164.762, 165.47, 162.502, 161.806, 161.14, 160.779, 164.834, 165.749, 164.855, 164.123, 166.515, 167.829, 168.756, 168.182, 170.023, 169.447, 170.372, 164.898, 163.997, 164.475, 165.616, 162.672, 163.559, 163.827, 163.588, 163.544, 163.406, 163.208, 161.781, 161.525, 164.028, 163.776, 
165.504, 160.853, 159.461, 158.56, 158.8, 159.271, 159.827, 159.276, 157.834, 156.917, 157.533, 156.605, 156.247, 156.009, 155.349, 155.599, 156.606, 154.786, 156.229, 155.9, 155.62, 156.232, 157.022, 154.684, 154.309, 154.899, 155.101, 152.794, 152.057, 150.642, 149.922, 150.249, 155.171, 155.738, 155.177, 155.232, 157.253, 154.635, 154.075, 154.579, 153.866, 152.546, 151.971, 155.807, 156.401, 156.418, 156.565, 157.822, 158.602, 159.116, 159.747, 156.266, 156.278, 157.644, 157.916, 155.21, 155.664, 154.921, 158.491, 159.853, 159.951, 161.041, 160.463, 159.741, 160.387, 160.161, 158.809, 158.795, 159.253, 159.66, 157.399, 159.206, 159.595, 160.695, 160.788, 158.377, 157.491, 156.444, 157.904, 161.526, 162.601, 163.805, 164.192, 163.057, 163.421, 161.953, 163.816, 164.376, 165.543, 166.361, 166.336, 166.459, 166.843, 165.937, 168.053, 167.092, 167.918, 167.147, 167.488, 168.324, 167.179, 166.107, 165.285, 164.978, 165.226, 163.958, 162.993, 164.235, 164.447, 164.118, 165.263, 165.107, 163.793, 163.538, 166.414, 167.569, 168.07, 168.547, 168.734, 169.346, 168.232, 167.969, 168.443, 167.619, 168.161, 168.513, 169.788, 169.867, 170.984, 166.31, 165.416, 165.287, 164.74, 164.032, 163.411, 164.12, 165.803, 165.758, 166.063, 165.476, 166.716, 166.352, 166.658, 166.975, 167.329, 166.793, 167.323, 168.839, 169.452, 169.335, 170.143, 169.898, 170.71, 170.589, 165.746, 165.156, 163.648, 163.1, 165.539, 167.03, 167.758, 167.725, 169.142, 169.115, 169.823, 171.206, 162.976, 161.525, 160.918, 159.72, 161.037, 161.673, 160.749, 160.138, 159.426, 159.236, 158.916, 161.748, 161.288, 162.226, 163.436, 161.282, 160.507, 159.183, 160.861, 158.754, 159.752, 161.671, 162.483, 161.748, 160.605, 162.843, 162.41, 161.826, 162.773, 163.999, 161.463, 160.914, 162.202, 163.004, 163.481, 162.684, 162.201, 161.884, 160.878, 163.152, 164.787, 165.355, 165.316, 164.9, 166.796, 166.893, 166.195, 166.633, 165.209, 165.745, 165.758, 165.296, 165.485, 167.173, 168.043, 164.702, 164.225, 164.034, 163.445, 162.876, 162.179, 162.84, 160.682, 164.509, 164.333, 164.66, 165.467, 165.206, 166.608, 167.241, 167.577, 168.562, 168.793, 167.536, 169.967, 168.71, 169.909, 164.035, 164.248, 165.185, 165.267, 162.904, 162.252, 163.141, 165.908, 166.843, 166.89, 166.71, 168.272, 169.204, 168.574, 167.129, 167.213, 168.553, 168.846, 166.079, 165.526, 164.963, 164.202, 169.359, 170.684, 170.656, 170.629, 171.616, 171.69, 170.571, 172.751, 170.94, 172.258, 170.671, 170.652, 171.885, 172.917, 170.741, 171.584, 170.965, 171.782, 172.914, 173.265, 172.704, 172.66, 171.654, 172.205, 174.202, 174.705, 173.65, 173.828, 175.52, 176.307, 176.953, 177.337, 172.561, 171.459, 170.648, 170.093, 170.535, 171.254, 171.985, 172.745, 171.755, 170.576, 169.837, 169.976, 168.981, 
170.315, 171.211, 171.367, 170.767, 170.23, 172.879, 173.306, 172.518, 170.845, 170.311, 168.835, 168.241, 171.12, 170.361, 171.471, 168.231, 166.835, 166.107, 165.931, 166.728, 167.259, 166.482, 168.536, 166.966, 169.032, 168.248, 165.677, 164.96, 163.925, 163.393, 164.326, 165.377, 165.789, 163.652, 162.67, 161.909, 162.504, 163.32, 164.057, 162.245, 160.589, 159.776, 159.684, 159.44, 158.357, 158.419, 157.475, 159.899, 159.831, 159.221, 159.574, 161.24, 162.208, 161.176, 158.283, 157.621, 157.725, 157.759, 157.797, 157.882, 156.887, 156.426, 159.275, 160.223, 159.792, 156.562, 155.617, 155.584, 155.242, 154.245, 153.02, 155.966, 155.955, 154.562, 153.954, 156.968, 158.349, 159.363, 158.915, 160.533, 160.285, 158.4, 161.149, 159.263, 160.62, 154.067, 152.714, 152.523, 153.173, 151.754, 150.461, 151.508, 151.627, 151.371, 150.815, 149.775, 150.313, 150.618, 150.906, 151.5, 151.058, 149.544, 148.957, 151.78, 148.914, 147.477, 146.735, 145.709, 147.135, 147.555, 146.542, 146.801, 145.369, 147.253, 146.649, 146.183, 146.882, 147.625, 147.072, 145.003, 144.379, 145.031, 144.597, 142.965, 143.17, 144.112, 146.067, 146.685, 147.477, 148.081, 147.55, 148.084, 148.69, 147.431, 148.099, 149.056, 148.816, 147.1, 147.221, 145.682, 150.127, 151.121, 150.542, 151.056, 151.978, 151.912, 150.442, 149.474, 148.871, 148.362, 147.825, 147.736, 148.54, 148.079, 149.038, 148.638, 146.718, 146.734, 145.355, 145.145, 144.392, 150.3, 151.286, 151.145, 151.019, 152.71, 153.013, 153.713, 154.431, 151.171, 151.037, 149.718, 149.488, 151.153, 150.571, 152.599, 148.847, 147.578, 147.791, 146.977, 146.489, 146.054, 144.928, 144.503, 143.451, 148.9, 149.22, 150.175, 150.896, 149.827, 148.869, 150.198, 150.176, 151.049, 152.524, 153.363, 150.698, 151.528, 151.495, 152.352, 152.261, 153.126, 153.074, 153.832, 152.839, 154.216, 154.588, 155.703, 153.636, 153.838, 155.223, 155.825, 155.728, 157.067, 157.371, 157.245, 158.084, 159.516, 159.729, 160.841, 158.664, 157.777, 158.1, 159.433, 159.869, 157.019, 156.006, 157.637, 
154.889, 160.077, 161.368, 161.507, 160.606, 162.478, 162.594, 162.96, 162.357, 163.091, 162.488, 162.856, 162.669, 162.95, 164.29, 165.04, 162.986, 162.906, 164.574, 165.816, 166.408, 165.743, 165.602, 165.348, 166.83, 166.53, 167.673, 168.31, 167.558, 167.098, 167.398, 166.72, 167.656, 167.794, 168.319, 169.257, 170.326, 170.037, 168.502, 171.573, 172.678, 172.666, 173.295, 172.599, 173.672, 172.636, 173.45, 171.979, 171.849, 173.017, 173.859, 170.545, 170.196, 171.167, 171.309, 171.844, 173.036, 174.028, 173.817, 173.219, 175.496, 176.327, 176.029, 174.33, 174.268, 174.079, 173.752, 173.108, 173.402, 172.111, 174.401, 174.297, 174.134, 172.714, 171.897, 174.368, 175.632, 172.407, 171.056, 170.052, 170.395, 171.232, 172.671, 173.356, 168.817, 167.761, 167.364, 167.29, 166.548, 166.721, 165.777, 166.472, 167.104, 166.706, 165.322, 164.894, 167.639, 168.985, 167.065, 168.919, 164.611, 163.283, 163.046, 162.639, 162.216, 160.855, 162.284, 163.334, 163.168, 161.717, 160.795, 163.564, 165.029, 165.977, 167.453, 168.414, 161.517, 160.181, 159.737, 160.452, 160.177, 158.681, 158.555, 158.082, 157.689, 156.506, 156.896, 156.336, 157.579, 158.688, 158.446, 157.618, 156.875, 159.759, 160.372, 161.628, 159.349, 157.747, 156.964, 155.564, 154.879, 157.534, 157.583, 158.423, 159.632, 157.868, 155.147, 153.841, 153.336, 152.431, 153.951, 152.642, 152.664, 153.63, 153.912, 153.528, 154.305, 155.504, 153.818, 155.248, 155.455, 155.072, 153.621, 154.234, 155.57, 
155.614, 153.276, 153.224, 154.538, 155.591, 154.48, 156.682, 158.029, 158.48, 159.628, 158.888, 158.226, 156.749, 157.603, 157.952, 157.759, 156.641, 157.047, 157.241, 156.711, 155.276, 154.7, 155.437, 153.386, 158.851, 158.731, 158.663, 158.432, 159.885, 159.962, 161.167, 158.885, 158.788, 158.39, 158.972, 160.106, 160.031, 161.202, 161.025, 162.061, 157.403, 156.966, 156.201, 156.133, 156.11, 156.1, 157.191, 155.014, 155.633, 154.864, 153.479, 152.79, 154.725, 155.981, 153.073, 151.771, 151.995, 153.092, 150.901, 151.703, 150.424, 151.027, 150.957, 151.057, 150.556, 149.366, 150.218, 150.307, 151.724, 152.209, 152.415, 151.433, 150.965, 150.509, 150.831, 152.054, 152.423, 151.543, 153.671, 151.905, 154.045, 153.166, 153.576, 149.776, 149.287, 149.895, 150.34, 147.783, 147.069, 145.592, 147.303, 149.906, 150.445, 151.961, 152.552, 150.014, 149.748, 150.583, 148.704, 152.576, 154.011, 154.368, 153.513, 154.429, 155.703, 155.637, 156.035, 155.964, 156.467, 157.463, 157.444, 156.685, 155.329, 155.165, 155.532, 155.037, 153.729, 153.574, 152.213, 151.139, 149.858, 156.373, 156.824, 156.04, 155.62, 158.303, 158.969, 158.954, 160.4, 155.858, 155.147, 156.067, 156.554, 153.918, 152.855, 151.635, 152.484, 156.307, 157.165, 156.337, 155.352, 158.309, 159.176, 159.147, 160.294, 156.73, 155.968, 156.843, 157.0, 154.862, 154.065, 157.409, 158.259, 157.453, 156.51, 159.313, 160.481, 159.811, 160.168, 157.834, 157.194, 158.34, 158.715, 156.214, 156.88, 155.031, 158.904, 160.046, 159.841, 158.71, 160.499, 160.965, 160.97, 160.102, 160.029, 162.391, 163.078, 162.459, 162.158, 162.273, 159.421, 158.589, 159.435, 158.991, 157.997, 159.046, 159.42, 160.662, 161.608, 162.629, 163.695, 162.369, 161.445, 163.116, 162.286, 163.12, 164.632, 165.176, 162.794, 165.326, 166.777, 167.436, 166.764, 166.973, 165.602, 164.712, 168.759, 169.474, 169.298, 168.818, 170.921, 170.774, 169.669, 169.689, 169.577, 168.727, 168.08, 170.961, 168.733, 167.953, 168.593, 169.818, 167.779, 167.079, 167.754, 168.238, 168.567, 169.218, 167.207, 165.892, 167.339, 168.127, 168.393, 168.129, 167.49, 167.523, 165.768, 168.624, 168.427, 167.771, 168.29, 169.748, 170.587, 169.49, 171.987, 166.627, 165.908, 166.194, 165.826, 164.42, 164.109, 163.692, 162.674, 166.844, 167.165, 166.523, 166.129, 168.687, 169.222, 169.361, 166.408, 165.83, 166.889, 167.454, 164.568, 164.213, 163.408, 167.169, 168.165, 167.544, 166.603, 169.299, 169.936, 168.089, 167.579, 168.399, 169.403, 167.961, 168.64, 167.606, 166.507, 169.268, 169.851, 170.344, 167.948, 167.011, 167.43, 168.616, 166.885, 166.111, 166.093, 164.695, 166.446, 166.707, 166.059, 164.834, 166.046, 166.436, 166.873, 166.334, 166.726, 167.537, 166.809, 166.632, 167.392, 169.061, 166.135, 166.432, 165.49, 164.363, 166.3, 164.889, 164.12, 164.07, 162.885, 162.827, 165.955, 165.254, 165.45, 164.697, 163.762, 163.529, 166.486, 166.814, 165.663, 164.911, 167.995, 167.642, 167.23, 165.541, 164.522, 165.215, 166.14, 163.941, 163.112, 162.164, 162.353, 164.773, 165.404, 165.159, 164.306, 164.969, 165.253, 165.724, 166.677, 165.931, 165.823, 167.128, 168.014, 165.56, 165.251, 166.776, 167.235, 168.434, 169.462, 169.177, 168.103, 169.33, 170.265, 169.362, 170.654, 171.704, 172.678, 173.779, 172.492, 172.711, 171.74, 172.271, 173.106, 173.671, 174.86, 172.271, 173.032, 172.801, 173.199, 172.439, 171.619, 172.858, 171.463, 173.688, 172.668]
})

fig = px.scatter_3d(df, x='X', y='Y', z='Z')

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),

    html.Div(children='''
        Dash: A web application framework for your data.
    '''),

    dcc.Graph(
        id='example-graph',
        figure=fig
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)