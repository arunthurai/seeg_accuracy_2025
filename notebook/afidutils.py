import numpy as np
import pandas as pd
import os
from glob import glob
import re
from itertools import cycle
from plotly.offline import plot
import plotly.graph_objs as go
import matplotlib.cm as cm
import random as rand
import nibabel as nib
import base64
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from io import BytesIO
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import pickle



# Dictionary for AFID labels
afids_labels = {
    1: 'AC', 2: 'PC', 3: 'ICS', 4: 'PMJ', 5: 'SIPF', 6: 'RSLMS',
    7: 'LSLMS', 8: 'RILMS', 9: 'LILMS', 10: 'CUL', 11: 'IMS', 12: 'RMB',
    13: 'LMB', 14: 'PG', 15: 'RLVAC', 16: 'LLVAC', 17: 'RLVPC', 18: 'LLVPC',
    19: 'GENU', 20: 'SPLE', 21: 'RALTH', 22: 'LALTH', 23: 'RSAMTH',
    24: 'LSAMTH', 25: 'RIAMTH', 26: 'LIAMTH', 27: 'RIGO', 28: 'LIGO',
    29: 'RVOH', 30: 'LVOH', 31: 'ROSF', 32: 'LOSF'
}

#hardcoding right and left hemisphere afids
right_afids = ['x_6','x_8','x_12','x_15','x_17','x_21','x_23','x_25','x_27','x_29','x_31','y_6','y_8','y_12','y_15','y_17','y_21','y_23','y_25','y_27','y_29','y_31', 'z_6','z_8','z_12','z_15','z_17','z_21','z_23','z_25','z_27','z_29','z_31']
left_afids = ['x_7','x_9','x_13','x_16','x_18','x_22','x_24','x_26','x_28','x_30','x_32','y_7','y_9','y_13','y_16','y_18','y_22','y_24','y_26','y_28','y_30','y_32', 'z_7','z_9','z_13','z_16','z_18','z_22','z_24','z_26','z_28','z_30','z_32']
combined_lables = ['AC','PC', 'ICS', 'PMJ','SIPF','SLMS','ILMS','CUL','IMS','MB','PG','LVAC','LVPC','GENU','SPLE','ALTH','SAMTH','IAMTH','IGO','VOH','OSF']
combined_lables = [element + axis for axis in ['x', 'y', 'z'] for element in combined_lables]
r = [6,8,12,15,17,21,23,25,27,29,31]
l = [7,9,13,16,18,22,24,26,28,30,32]

exclude_afids_labels = [10,14,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
exclude_afids = ['CULx','PGx','GENUx', 'SPLEx', 'ALTHx', 'SAMTHx', 'IAMTHx', 'IGOx','VOHx','OSFx',
                 'CULy','PGy','GENUy', 'SPLEy', 'ALTHy', 'SAMTHy', 'IAMTHy', 'IGOy','VOHy','OSFy',
                 'CULz','PGz','GENUz', 'SPLEz', 'ALTHz', 'SAMTHz', 'IAMTHz', 'IGOz','VOHz','OSFz'
                 ]

def fcsvtodf(fcsv_path):
    """
    Convert a .fcsv file (assumes RAS coordinate system) to a ML-friendly dataframe and return the cleaned xyz coordinates.
    
    Parameters:
    - fcsv_path: str, path to the .fcsv file
    
    Returns:
    - df_xyz_clean: pandas.DataFrame with cleaned x, y, z coordinates
    - num_points: int, number of fiducial points
    """

    # Extract the subject ID from the file path (naming is in bids-like)
    subject_id = re.search(r'(sub-\w+)', fcsv_path).group(1)

    # Read in .fcsv file, skip header
    df_raw = pd.read_table(fcsv_path,sep=',',header=2)

    #Extract the x, y, z coordiantes and store them in data science friendly format (i.e., features in cols and subject in rows)
    df_xyz = df_raw[['x','y','z']].melt().transpose()

    #Use number of row in fcsv to make number points
    colnames = [f'{axis}_{i % int(df_raw.shape[0]) + 1}' for axis in ['x', 'y', 'z'] for i in range(int(df_raw.shape[0]))]
    
    #Reassign features to be descriptive of coordinate
    df_xyz.columns = colnames

    #clean dataframe and pin correct subject name
    df_xyz_clean = df_xyz.drop('variable', axis= 0)
    df_xyz_clean = df_xyz_clean.rename(index={'value': subject_id})
    df_xyz_clean = df_xyz_clean.astype(float)

    return df_xyz_clean, df_raw.shape[0]


def get_fiducial_index(fid):
    """
    Retrieve the index corresponding to the fiducial name or integer.
    
    Parameters:
    - fid: str or int, fiducial identifier (name or index)
    
    Returns:
    - int, corresponding fiducial index
    """

    if isinstance(fid, str):
        for idx, name in afids_labels.items():
            if name == fid:
                return idx
    elif isinstance(fid, int):
        return fid
    raise ValueError("Invalid fiducial identifier.")


def compute_distance(fcsv_path, fid1, fid2):
    """
    Compute the Euclidean distance between two fiducials.

    Parameters:
    - fcsv_path: str, path to the .fcsv file
    - fid1, fid2: str or int, fiducial identifiers

    Returns:
    - xyz_diff: numpy.array, difference in x, y, z coordinates
    - distance: float, Euclidean distance between fiducials
    """

    # Retrieve indices of the fiducials
    index1, index2 = get_fiducial_index(fid1), get_fiducial_index(fid2)

    # Load dataframe from the fcsv file
    df = fcsvtodf(fcsv_path)[0]

    # Extract x, y, z coordinates into numpy arrays
    coords1 = df[[f'x_{index1}', f'y_{index1}', f'z_{index1}']].to_numpy()
    coords2 = df[[f'x_{index2}', f'y_{index2}', f'z_{index2}']].to_numpy()

    # Compute the difference as a numpy array
    xyz_diff = coords1 - coords2

    # Compute the Euclidean distance
    distance = np.linalg.norm(xyz_diff)

    return xyz_diff.flatten(), distance


def compute_average(fcsv_path, fid1, fid2):
    """
    Compute the average position between two fiducials.
    
    Parameters:
    - fcsv_path: str, path to the .fcsv file
    - fid1, fid2: str or int, fiducial identifiers
    
    Returns:
    - xyz_average: numpy.array, average coordinates (x, y, z) between fiducials
    """

    # Retrieve indices of the fiducials
    index1, index2 = get_fiducial_index(fid1), get_fiducial_index(fid2)

    # Load dataframe from the fcsv file
    df = fcsvtodf(fcsv_path)[0]

    # Extract x, y, z coordinates into numpy arrays
    coords1 = df[[f'x_{index1}', f'y_{index1}', f'z_{index1}']].to_numpy()
    coords2 = df[[f'x_{index2}', f'y_{index2}', f'z_{index2}']].to_numpy()

    # Compute the average as a numpy array
    xyz_average = (coords1 + coords2)/2

    return xyz_average.flatten()


def acpcmatrix(fcsv_path, midline, center_on_mcp = True, write_matrix = False, transform_file_name = None):
    """
    Computes a 4x4 transformation matrix aligning with the AC-PC axis.
    
    Parameters:
    - fcsv_path: str, path to the .fcsv file
    - midline: str or int, one of the 10 midline AFID points.
    - center_on_mcp: bool, If True, adds translation element to the ACPC matrix (centering on MCP).
    
    Returns:
    - matrix: np.ndarray, A 4x4 affine transformation matrix.
    """

    # A-P axis
    acpc = compute_distance(fcsv_path, 'AC', 'PC')  # vector from PC to AC
    yAxis = acpc[0] / acpc[1]  # unit vector defining anterior and posterior axis

    # R-L axis
    lataxis = compute_distance(fcsv_path, midline, 'AC')[0] # vector from AC to midline point 
    xAxis = np.cross(yAxis, lataxis) # vector defining left and right axis
    xAxis /= np.linalg.norm(xAxis) # unit vector defining left and right axis

    # S-I axis
    zAxis = np.cross(xAxis, yAxis)
    zAxis /= np.linalg.norm(zAxis)

    # Rotation matrix
    rotation = np.vstack([xAxis, yAxis, zAxis])
    translation = np.array([0, 0, 0])

    # Build 4x4 matrix
    matrix = np.eye(4)
    matrix[:3, :3] = rotation
    matrix[:3, 3] = translation

    # Orientation correction for matrix is midpoint placed below ACPC (i.e., at PMJ)
    if matrix[0][0] < 0:
        matrix = matrix * np.array([[-1,-1,-1,0],
                                    [1,1,1,0],
                                    [-1,-1,-1,0],
                                    [0,0,0,1]])

    if center_on_mcp: #need to compute MCP AFTER rotation is applied 
        mcp = compute_average(fcsv_path,'AC','PC') #MCP in native 
        matrix[:3, 3] = -np.dot(matrix[:3, :3], mcp) # MCP after rotation; negative because we set MCP to (0,0,0)

    if write_matrix: 
        generate_slicer_file(matrix, transform_file_name)

    return matrix

def extract_sub_metadata(fcsv_path):
    """
    Extracts subject id, description, and dataset following bids-like structure.
    
    Parameters:
    - fcsv_path: str, path to the .fcsv file
    
    Returns:
    - sub_id: str, subject ID 
    - description: str, description label
    - dataset: str, dataset
    """    

    filename = os.path.basename(fcsv_path)

    sub_id_match = re.search(r'^sub-[A-Za-z0-9]+', filename)
    sub_id = sub_id_match.group() if sub_id_match else 'unknown'

    description_match = re.search(r'desc-([A-Za-z0-9\-]+)', filename)
    description = description_match.group(1) if description_match else 'unknown'
    
    parts = fcsv_path.split(os.sep)

    dataset = "unknown_dataset"

    for i, part in enumerate(parts):
        if part == 'anat':
            if i >= 3:
                dataset = parts[i - 3]
            break  # exit after finding the first match


    match = re.search(r'desc-[a-zA-Z0-9_]+_([a-zA-Z0-9]+)\.fcsv', parts[-1])
    if match:
        coord_type = match.group(1) 
    else: 
        coord_type = "unknown_coordinatetype"

    return sub_id, description, dataset, coord_type

    
def transform_afids(fcsv_path, midline_fcsv=None, output_dir=None, midpoint='PMJ', save_matrix = False):
    """
    Computes and applies an AC-PC transformation to AFID files and saves an fcsv with transfromed coordinates.

    Parameters:
    - fcsv_path: str, path to the .fcsv file for coordinates to be transformed
    - midline_fcsv: str, path to the .fcsv file which has the midline coordinates
    - output_dir: str, path of the directory to store transformed .fcsv file (if not specified, no output fcsv will be written)
    - midpoint: str or int, any midline AFID point
   
    Returns:
    - tcoords: np.ndarray, transformed coordinates
    """
    if midline_fcsv:
        # Compute the 4x4 AC-PC transformation matrix
        xfm_txt = acpcmatrix(midline_fcsv, midpoint)
    else: 
        xfm_txt = acpcmatrix(fcsv_path, midpoint)
        
    # Read header lines separately
    with open(fcsv_path, 'r') as file:
        header_lines = [next(file) for _ in range(3)]

    # Read coordinates from the file
    fcsv_df = pd.read_table(fcsv_path, sep=",", header=2)
    
    # Copy coordinates and apply transformation
    pretfm = fcsv_df.loc[:, ["x", "y", "z"]].copy()
    pretfm["coord"] = 1  # Add a fourth column for homogeneous coordinates
    coords = pretfm.to_numpy()
    thcoords = xfm_txt @ coords.T
    tcoords = thcoords.T[:, :3]  # Remove homogeneous coordinates
    
    # Substitute new coordinates
    if tcoords.shape == (len(fcsv_df), 3):
        fcsv_df.loc[:, "x"] = tcoords[:, 0]
        fcsv_df.loc[:, "y"] = tcoords[:, 1]
        fcsv_df.loc[:, "z"] = tcoords[:, 2]
    else:
        raise ValueError("New coordinates do not match the number of rows in the original fcsv.")
    
    if output_dir is not None:
        # Save the transformed coordinates to a new fcsv file
        sub_desc_dataset = extract_sub_metadata(fcsv_path)
        output_filename = f'{sub_desc_dataset[0]}_space-MCP{midpoint}_desc-{sub_desc_dataset[1]}_{sub_desc_dataset[3]}.fcsv'
        output_bids = f"{sub_desc_dataset[2]}/{sub_desc_dataset[3]}_mcp/{sub_desc_dataset[0]}/anat"
        full_out_dir = os.path.join(output_dir, output_bids)
        if not os.path.exists(full_out_dir):
            os.makedirs(full_out_dir)

        output_path = os.path.join(full_out_dir, output_filename)
        
        # Write the header and updated DataFrame to the new fcsv file with specific formatting
        with open(output_path, 'w') as file:
            file.writelines(header_lines)
            fcsv_df.to_csv(file, sep=",", index=False, header=False, float_format="%.5f")
        if save_matrix:
            # Reconstruct path with approriate extension 
            tfm_path = os.path.splitext(output_path)[0] + '.txt'
            generate_slicer_file(xfm_txt,tfm_path)
    
        print(f"Saved transformed file: {output_filename}")
    
    return tcoords


def coords_to_fcsv(coords, fcsv_template, fcsv_output):
    """
    Updates an fcsv template file with new coordinates and writes the updated version to a new file.
    Ensures the directory for the output file exists before writing.

    Parameters:
    - coords: ndarray, array of shape (n, 3) containing the coordinates for each fiducial.
    - fcsv_template: str, path to a template fcsv file that can be used as a reference.
    - fcsv_output: str, path where the output fcsv file with updated coordinates will be written.
    """

    # Ensure output directory exists
    output_dir = os.path.dirname(fcsv_output)
    os.makedirs(output_dir, exist_ok=True)

    # Read in the fcsv template
    with open(fcsv_template, "r") as f:
        fcsv = [line.strip() for line in f]

    # Loop over each fiducial in the coordinates
    for fid in range(1, coords.shape[0] + 1):
        # Update fcsv, skipping header lines (first 3 lines assumed to be header)
        line_idx = fid + 2  # Line to update with new coordinates
        centroid_idx = fid - 1  # Index of the current point in coords

        # Replace the placeholders for x, y, and z coordinates
        fcsv[line_idx] = fcsv[line_idx].replace(
            f"afid{fid}_x", str(coords[centroid_idx][0])
        ).replace(
            f"afid{fid}_y", str(coords[centroid_idx][1])
        ).replace(
            f"afid{fid}_z", str(coords[centroid_idx][2])
        )

    # Write the updated fcsv to the output file
    with open(fcsv_output, "w") as f:
        f.write("\n".join(fcsv))


    
def generate_slicer_file(matrix, filename):
    """
    Generate a .txt transformation file for 3D Slicer from a 4x4 matrix.

    Parameters: 
    - matrix: np.ndarray, 4x4 transformation matrix
    - output_path: str, path to store .txt file
    """
    D = np.array([
    [-1,  0,  0, 0],
    [ 0, -1,  0, 0],
    [ 0,  0,  1, 0],
    [ 0,  0,  0, 1]
    ])

    ras_inmatrix = np.linalg.inv(matrix)

    lps_inmatrix = D @ ras_inmatrix @ D

    # Extract rotation/scale and translation components
    rotation_scale = lps_inmatrix[0:3, 0:3].flatten()
    translation = lps_inmatrix[0:3, 3]

    # Format the content of the .tfm file
    tfm_content = "#Insight Transform File V1.0\n"
    tfm_content += "#Transform 0\n"
    tfm_content += "Transform: AffineTransform_double_3_3\n"
    tfm_content += "Parameters: " + " ".join(map(str, rotation_scale)) + " " + " ".join(map(str, translation)) + "\n"
    tfm_content += "FixedParameters: 0 0 0\n"

    # Write the content to the specified file
    with open(filename, 'w') as file:
        file.write(tfm_content)

    print(f"{filename} has been generated.")


def plot_multiple_files(file_list, out_file, rotation=False, trace_lines=True, save_gif=True, gif_filename='rotation.gif'):
    traces = []
    all_afids = []
    all_points = pd.DataFrame()

    # Gather all AFIDs from all files for consistent coloring
    for file_path in file_list:
        data = pd.read_csv(file_path, skiprows=2)
        subject_id, dataset_name = extract_sub_metadata(file_path)[0], extract_sub_metadata(file_path)[2]

        if '_stn.fcsv' in file_path:
            afid_type = '_stn'
            data['file'] = f"{subject_id}{afid_type}"
            data['dataset'] = f"{dataset_name}{afid_type}"
            data['label'] = data['desc']
            all_afids.extend(data['desc'].unique())
        else:
            data['file'] = subject_id
            data['dataset'] = dataset_name
            all_afids.extend(data['label'].unique())

        all_points = pd.concat([all_points, data], ignore_index=True)

    # Get unique AFIDs and assign colors
    unique_afids = list(set(all_afids))
    rand.shuffle(unique_afids)
    color_map = cm.get_cmap('tab20', len(unique_afids))  # Use tab20 colormap for a wide range
    afid_colors = {afid: color_map(i) for i, afid in enumerate(unique_afids)}

    # Convert colors to RGBA format
    afid_colors = {afid: f'rgba({int(r * 255)}, {int(g * 255)}, {int(b * 255)}, {a})'
                   for afid, (r, g, b, a) in afid_colors.items()}

    # Calculate centroid for each AFID and keep the descriptions
    centroids = all_points.groupby('label').agg({'x': 'mean', 'y': 'mean', 'z': 'mean'}).reset_index()

    # Get the unique dataset names
    datasets = all_points['dataset'].unique()

    # Create a dummy trace for each dataset to allow toggling all subjects in the dataset
    for dataset in datasets:
        dataset_trace = go.Scatter3d(
            x=[None], y=[None], z=[None],  
            mode='markers',
            marker=dict(size=1),
            showlegend=True,
            name=f'{dataset}',
            legendgroup=dataset,
            hoverinfo='none'
        )
        traces.append(dataset_trace)

    # Create traces for original points, grouped by dataset
    for file_path in file_list:
        subject_id = extract_sub_metadata(file_path)[0]
        afid_type = '_stn' if '_stn.fcsv' in file_path else ''
        data = all_points[all_points['file'] == f"{subject_id}{afid_type}"]

        x_coords = data['x']
        y_coords = data['y']
        z_coords = data['z']
        afids = data['label']
        dataset_name = data['dataset'].iloc[0]

        # Get colors for each point based on its AFID
        colors = [afid_colors[afid] for afid in afids]

        # Create a scatter trace for each subject, grouped by dataset
        trace = go.Scatter3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            mode='markers',
            marker=dict(
                size=6,  
                color=colors,
                opacity=1, 
                line=dict(width=0.5, color='black') 
            ),
            name=subject_id,
            legendgroup=dataset_name,
            legendgrouptitle=dict(text=dataset_name),
            showlegend=True
        )
        traces.append(trace)

    # Add centroid markers with the AFID descriptions as labels
    centroid_trace = go.Scatter3d(
        x=centroids['x'],
        y=centroids['y'],
        z=centroids['z'],
        mode='markers+text',
        text=centroids['label'],  # Label with the AFID descriptions
        textposition='top center',
        marker=dict(
            size=8,  # Larger marker for centroids
            color='black',
            opacity=1,
            symbol='diamond',
            line=dict(width=2)
        ),
        name='Centroid',
        showlegend=True
    )
    traces.append(centroid_trace)

    if trace_lines:
        # Add lines from all centroids to the STN centroids
        stn_centroids = centroids[centroids['label'].str.contains('_stn')]

        # Iterate through each STN centroid
        for _, stn_row in stn_centroids.iterrows():
            # Draw lines from all centroids to the current STN centroid
            for _, centroid_row in centroids.iterrows():
                line_trace = go.Scatter3d(
                    x=[centroid_row['x'], stn_row['x']],
                    y=[centroid_row['y'], stn_row['y']],
                    z=[centroid_row['z'], stn_row['z']],
                    mode='lines',
                    line=dict(color='grey', width=2),  # Use a blue line to distinguish these connections
                    showlegend=True  # No need to add to the legend
                )
                traces.append(line_trace)

    # Set layout with more visual enhancements
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X Coordinate', backgroundcolor='rgb(240, 240, 240)', showgrid=True, gridcolor='white'),
            yaxis=dict(title='Y Coordinate', backgroundcolor='rgb(240, 240, 240)', showgrid=True, gridcolor='white'),
            zaxis=dict(title='Z Coordinate', backgroundcolor='rgb(240, 240, 240)', showgrid=True, gridcolor='white'),
            aspectmode='data'
        ),
        title="3D Interactive Plot of Fiducial Points with Centroids",
        margin=dict(l=0, r=0, b=0, t=30),
        legend=dict(title="Legend", itemsizing='constant', groupclick="toggleitem"),
        scene_camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.5)
        ),
        hovermode='closest'
    )

    # Create the figure
    fig = go.Figure(data=traces, layout=layout)

    if rotation:
        # Prepare animation frames for a smoother 360-degree rotation
        frames = [go.Frame(layout=dict(scene_camera=dict(eye=dict(
            x=1.5 * np.cos(np.radians(angle)),
            y=1.5 * np.sin(np.radians(angle)),
            z=1.5))))
            for angle in range(0, 360, 10)]  # Reduced frame step for smoother animation

        # Add the play button for the rotation animation
        fig.update_layout(
            updatemenus=[dict(type='buttons', showactive=False,
                              buttons=[dict(label='Play',
                                            method='animate',
                                            args=[None, dict(frame=dict(duration=100, redraw=True), fromcurrent=True)])])],
        )

        # Attach frames to the figure
        fig.frames = frames

    # Display the interactive plot
    plot(fig, filename=out_file)

    if save_gif:
        # Save the rotating plot as a gif
        colors_per_afid = {label: (np.random.rand(), np.random.rand(), np.random.rand(), 1) for label in unique_afids}
        def update(angle):
            ax.view_init(elev=15, azim=angle)
            ax.clear()

            # Plot the points with different colors for each label
            for label in unique_afids:
                label_points = all_points[all_points['label'] == label]
                ax.scatter(label_points['x'], label_points['y'], label_points['z'],
                           s=10, alpha=0.8, linewidths=0.7, edgecolors='k',
                           color=colors_per_afid[label], label=label)

            # Plot centroids with improved visuals
            ax.scatter(centroids['x'], centroids['y'], centroids['z'],
                       s=10, c='black', alpha=1, marker='D', edgecolors='k', label='Centroid')

            # Set axis labels for better context
            ax.set_xlabel('X Coordinate')
            ax.set_ylabel('Y Coordinate')
            ax.set_zlabel('Z Coordinate')

        fig_mpl = plt.figure(figsize=(10, 8))
        ax = fig_mpl.add_subplot(111, projection='3d')
        ani = FuncAnimation(fig_mpl, update, frames=np.arange(0, 360, 2), interval=50)
        ani.save(gif_filename, writer='pillow')

    return out_file

def save_single_slice_in_memory(data, x, y, z, offset, zoom_radius, show_crosshairs):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (slice_data, coord, title) in enumerate(
        zip(
            [data[x + offset, :, :], data[:, y + offset, :], data[:, :, z + offset]],
            [(y, z), (x, z), (x, y)],
            ["Sagittal", "Coronal", "Axial"],
        )
    ):
        axes[i].imshow(slice_data.T, origin='lower', cmap='gray')
        if offset == 0 and show_crosshairs:
            axes[i].axhline(y=coord[1], color='r', lw=1)
            axes[i].axvline(x=coord[0], color='r', lw=1)
        axes[i].set_xlim(coord[0] - zoom_radius, coord[0] + zoom_radius)
        axes[i].set_ylim(coord[1] - zoom_radius, coord[1] + zoom_radius)
        axes[i].set_title(title)
        axes[i].axis('off')

    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    plt.close(fig)
    buffer.seek(0)

    return base64.b64encode(buffer.read()).decode()

def save_mri_slices_as_images(data, x, y, z, jitter, zoom_radius, show_crosshairs=True):
    return Parallel(n_jobs=-1)(
        delayed(save_single_slice_in_memory)(data, x, y, z, offset, zoom_radius, show_crosshairs)
        for offset in range(-jitter, jitter + 1)
    )

def extract_coordinates_from_fcsv(file_path, label_description):
    df = pd.read_csv(file_path, comment='#', header=None)
    row = df[df[11] == label_description]
    return tuple(row.iloc[0, 1:4]) if not row.empty else None

def generate_html_with_keypress(subject_images, reference_images, output_html="mri_viewer.html"):
    """Generate an interactive HTML viewer with sticky instructions."""
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>MRI Viewer</title>
        <style>
            body { display: flex; font-family: Arial, sans-serif; margin: 0; padding: 0; }
            .instructions {
                width: 20%; padding: 20px; background-color: #f4f4f4;
                border-right: 2px solid #ddd;
                position: sticky; top: 0; height: 100vh; overflow-y: auto;
                box-sizing: border-box;
            }
            .instructions h2 { margin-top: 0; }
            .viewer { flex: 1; padding: 20px; text-align: center; }
            .slider { width: 80%; margin: 20px auto; }
            .image { display: block; margin: 0 auto; max-width: 80%; transition: opacity 0.2s ease-in-out; }
        </style>
    </head>
    <body>
        <div class="instructions">
            <h2>Tips & Tricks</h2>
            <ul>
                <li>Use the <strong>slider</strong> to navigate through MRI. Keyboard can also be used after clicking on slider.</li>
                <li>Check the <strong>Red</strong> crosshair which represents the placement of a given landmark.</li>
                <li>Press <strong>R</strong> to toggle between the subject scan and the protocol-defined placement (if one exists).</li>
                <li>Use the <strong>TAB</strong> key to navigate to next slider & MRI.</li>          
            </ul>
        </div>
        <div class="viewer">
    """
    
    for label, images in subject_images.items():
        num_slices = len(images)
        has_reference = reference_images is not None and label in reference_images

        html_content += f"""
        <div class="container">
            <h2>Landmark: {label}</h2>
            <img id="{label}_image" class="image" src="data:image/png;base64,{images[0]}" alt="MRI Slice">
            <input type="range" min="0" max="{num_slices - 1}" value="0" class="slider" id="{label}_slider">
        </div>
        """
    
    html_content += """
        </div>
        <script>
            const subjects = {};
    """
    
    for label, images in subject_images.items():
        has_reference = reference_images is not None and label in reference_images
        html_content += f"""
            subjects["{label}"] = {{
                targetImages: {images},
                referenceImages: {reference_images[label] if has_reference else "null"},
                showingReference: false
            }};
        """
    
    html_content += """
            document.addEventListener('DOMContentLoaded', () => {
                for (const [label, data] of Object.entries(subjects)) {
                    const slider = document.getElementById(`${label}_slider`);
                    const image = document.getElementById(`${label}_image`);

                    slider.addEventListener('input', () => updateImage(label));
                    document.addEventListener('keydown', (event) => {
                        if (event.key.toLowerCase() === 'r' && data.referenceImages) {
                            data.showingReference = !data.showingReference;
                            updateImage(label);
                        }
                    });

                    function updateImage(label) {
                        const sliceIndex = document.getElementById(`${label}_slider`).value;
                        const imageArray = subjects[label].showingReference
                            ? subjects[label].referenceImages
                            : subjects[label].targetImages;
                        document.getElementById(`${label}_image`).src = 'data:image/png;base64,' + imageArray[sliceIndex];
                    }
                }
            });
        </script>
    </body>
    </html>
    """
    
    with open(output_html, "w") as f:
        f.write(html_content)
        
def generate_interactive_mri_html(nii_path, fcsv_path, labels, ref_nii_path=None, ref_fcsv_path=None, jitter=2, zoom_radius=20, out_file_prefix="mri_viewer"): 
    """
    Generates an interactive HTML viewer for MRI slices based on fiducial coordinates for a single subject.

    Parameters:
    - nii_path: str
        Full path to the subject's NIfTI (.nii.gz) file.
    - fcsv_path: str
        Full path to the subject's FCSV (.fcsv) file.
    - labels: list of int
        Landmark indices to extract for the subject.
    - ref_nii_path: str or None
        Full path to reference NIfTI file (optional).
    - ref_fcsv_path: str or None
        Full path to reference FCSV file (optional).
    - jitter: int, optional (default=2)
        The number of pixels to expand around the coordinate in each slice for visualization.
    - zoom_radius: int, optional (default=20)
        The radius (in pixels) around the coordinate to extract for display.
    - out_file_prefix: str, optional (default="mri_viewer")
        The prefix for the output HTML file name.
    
    Notes:
    - Extracts specified fiducial coordinates from the subject's .fcsv file.
    - Maps the coordinates from world space to voxel space using the NIfTI affine transformation.
    - MRI slices centered around these coordinates are saved as images.
    - If a reference MRI and coordinate file are provided, the same process is applied.
    - Generates an interactive HTML viewer allowing keypress navigation.
    
    Returns:
    - None (outputs an HTML file with the interactive MRI viewer for the subject).
    """
    subject_key = os.path.basename(nii_path).replace(".nii.gz", "")
    target_img = nib.load(nii_path, mmap=True)
    affine_inv = np.linalg.inv(target_img.affine)
    target_data = target_img.get_fdata(dtype=np.float32)

    target_images = {}  # Dictionary to store images per landmark
    
    for label in labels:
        target_coord = extract_coordinates_from_fcsv(fcsv_path, label)
        if not target_coord:
            print(f"Coordinates for label '{label}' not found in subject '{subject_key}'.")
            continue
        
        target_voxel = np.round(affine_inv.dot(target_coord + (1,))).astype(int)
        target_images[label] = save_mri_slices_as_images(target_data, *target_voxel[:3], jitter, zoom_radius)
        
    reference_images = None
    if ref_nii_path and ref_fcsv_path:
        ref_img = nib.load(ref_nii_path, mmap=True)
        ref_data = ref_img.get_fdata(dtype=np.float32)
        reference_images = {}
        
        for label in labels:
            ref_coord = extract_coordinates_from_fcsv(ref_fcsv_path, label)
            if ref_coord:
                ref_voxel = np.round(np.linalg.inv(ref_img.affine).dot(ref_coord + (1,))).astype(int)
                reference_images[label] = save_mri_slices_as_images(ref_data, *ref_voxel[:3], jitter, zoom_radius)
    os.makedirs(out_file_prefix, exist_ok=True)

    out_file = f"{out_file_prefix}{subject_key}.html"
    
    print(f"Generated: {out_file}")

    # Pass the correct image dictionaries
    generate_html_with_keypress(target_images, reference_images, out_file)

def flip_and_concatenate_hemispheres(df_afids):
    """
    Flips hemispheres by removing left/right AFIDs, relabels columns, and concatenates 
    mirrored and non-mirrored data.

    Parameters:
        df_afids (pd.DataFrame): Original DataFrame of AFID coordinates.
        afidutils (module or object): Utility with attributes:
            - left_afids (list): Column names for left hemisphere AFIDs.
            - right_afids (list): Column names for right hemisphere AFIDs.
            - combined_lables (list): Unified column names after hemisphere flip.
            - exclude_afids (list): Columns to exclude from final DataFrame.

    Returns:
        pd.DataFrame: Combined and cleaned DataFrame with flipped hemispheres.
    """
    df_afids_mcp_l = df_afids.copy()
    df_afids_mcp_l[[col for col in df_afids.columns if 'x' in col]] *= -1
    
    # Drop left hemisphere from one copy, right from the other
    df_afids_mcp = df_afids.drop(left_afids, axis=1)
    df_afids_mcp_l = df_afids_mcp_l.drop(right_afids, axis=1)

    # Rename columns to have consistent structure
    df_afids_mcp.columns = combined_lables
    df_afids_mcp_l.columns = combined_lables

    # Concatenate mirrored and original data
    concatenated_df = pd.concat([df_afids_mcp, df_afids_mcp_l])

    # Drop excluded AFIDs
    concatenated_df.drop(exclude_afids, axis=1, inplace=True)

    return concatenated_df


def predict_stn_coordinates(
    subject_afid_path,
    model_path
):
    """
    Runs the trained AFIDs -> STN model on a single subject and returns predicted coordinates.

    Parameters
    ----------
    subject_afid_path : str
        Path to the subject's .fcsv AFID file in ACPC + mcp centered space.
    model_path : str
        Path to the trained model pickle file.
    afidutils : module
        AFID utility module, must include:
        - fcsvtodf()
        - left_afids, right_afids, combined_lables, exclude_afids

    Returns
    -------
    y_pred : np.ndarray
        Predicted STN coordinates, shape (82, 3).
    """

    # Load subject AFIDs
    df = fcsvtodf(subject_afid_path)[0]

    # values between -0.01 and 0.01 are set to zero
    df.map(make_zero) 

    df_proc = flip_and_concatenate_hemispheres(df)

    # Load model
    with open(model_path, 'rb') as file:
        model_dict = pickle.load(file)

    scaler = model_dict['standard_scaler']
    pca = model_dict['pca']
    ridge_models = [model_dict['x'], model_dict['y'], model_dict['z']]

    # Transform
    X = scaler.transform(df_proc.values)
    X = pca.transform(X)

    # Predict STN
    y_pred = np.column_stack([ridge.predict(X) for ridge in ridge_models])

    return y_pred

def make_zero(num, threshold=0.01):
    if abs(num) < threshold:
        return 0
    else:
        return num

def remove_outliers(df, cols, threshold=3.5):
    """
    Remove rows from a DataFrame based on a modified Z-score for a given column(s).
    """
    for col in cols:
        median = np.median(df[col])
        mad = np.median(np.abs(df[col] - median))
        if mad == 0:
            mad = 1.0
        z_scores = 0.6745 * (df[col] - median) / mad
        df = df.loc[abs(z_scores) < threshold]
    return df

def create_list(n):
    numbers = list(range(1, n + 1))  # Create a list of numbers from 1 to n
    repeated_numbers = cycle(numbers)  # Create a cycle iterator for the numbers
    result = [next(repeated_numbers) for _ in range(n * 2)]  # Generate n*2 elements by cycling through the numbers
    return result