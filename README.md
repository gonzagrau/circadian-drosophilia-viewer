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

**Prerequisites:** `python >= 3.11.9`, `git lfs`

First, if you don't have git lfs set up, please do so by running:

*On Linux*:

```
sudo apt-get install git-lfs
git lfs install
```

*On MacOS*:
```
brew install git-lfs
git lfs install
```

*On Windows*:
Git LFS should be downloaded by default, so we simply activate it by running
```
git lfs install
```

**To download this demo**, simply clone the repository

Using HTTPS:
```
git clone https://github.com/gonzagrau/circadian-drosophilia-viewer
```

Then, install the required packages. We reccomend using a virtual environment. For instance, with conda:
```
conda create -n circdrosvenv python=3.11.1
conda activate circdrosvenv
```

Using pip, get the requirements
```
pip install -r requirements.txt
```

Finally, to run the streamlit app, use:
```
streamlit run streamlit_main.py
```




