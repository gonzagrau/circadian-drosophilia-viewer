# Circadian Drosophilia Viewer
## An interactive visualizer for Rosbash's _transcriptomic taxonomy of Drosophila circadian neurons around the clock (2021)_

### Step 1: select the genes you want to visualize
![image](https://github.com/user-attachments/assets/547b204f-b112-4ddd-a912-d7a440e4cc7a)

### Step 2: download the normalized, annotated data as a csv
![image](https://github.com/user-attachments/assets/0de6ca2e-b159-4d3e-a12e-9d6946b85aa6)

### Step 3: Design your plots, and visualize!
![image](https://github.com/user-attachments/assets/64364db2-2401-4250-b7ff-4ec31cfc8a9c)

![image](https://github.com/user-attachments/assets/7d3de4bf-3a92-49ec-b26c-c4351dbbf9e2)

## Installation guide

Prerequisites: `python >= 3.11.9`

We reccommend creating a virtual environment for this project. For instance, using conda:
```
conda create -n circdrosvenv python=3.11.1
conda activate circdrosvenv
```

To download this demo, simply clone the repository

Using SSH:
```
git clone git@github.com:gonzagrau/circadian-drosophilia-viewer.git
```

Then, install the required packages:

```
pip install -r requirements.txt
```

Finally, to run the streamlit app, use:
```
streamlit run streamlit_main.py
```




