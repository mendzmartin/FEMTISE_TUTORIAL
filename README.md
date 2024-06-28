# Setting up an Organised Folder Structure for Research Projects
Here is a brief explanation of the folder hierarchy. Note that every folder is numbered in order to preserve a fixed ordering.
+ `010_project_managment`: obviously enough, this is the folder where you keep all your files related to managing and planning your research project. For example: background, proposals and plans, funding applications, budget, reports.
+ `020_ethics_governance`: in cognitive science, more so than perhaps in many other sciences, meeting legislative requirements and properly managing ethical and privacy issues is essential. This is why this folder structure contains a separate folder for documents such as: applications for ethical approval, insurance, participant information sheets, and consent forms.
+ `030_experimentation`: this folder is the heart of your research project and contains all files related to the design and implementation of your experiment (or simulation), as well as data collection, processing, and analysis. The experiment folder contains several subfolders to help you organise your files:
  + `010_input`: this is the place to store any files and documents used in the experimental task itself. This could be anything from picture and word stimuli, audio recordings an task instructions.
  + `020_src`: is the folder where you keep main software code used to run experiment, also any other run script or auxiliary code.
  + `030_test`: here you need to save input, output and software code to testing main software code, also to keep testing analysis output or calculations.
  + `040_output`: this folder contains all the raw data files collected during the experimental session. Only original data files should be kept here.
  + `050_analysis`: here you can keep all “intermediary” files: files pre-processed or transformed for further statistical analyses, including any analysis scripts or notebooks. Also this folder is useful to store all the analysis “output”, also graphical output such as figures and diagrams. Having these files in one place will make the experiment write up easier, for example for reports, conference presentations, or publications.
+ `040_dissemination`: the last but not the least important folder is the dissemination folder. As the name says, this is a central location to keep organised any documents related to disseminating your research findings to the scientific and wider community. This is where your journal articles will reside, but also any conference papers or posters which report the experimental findings. Not to be dismissed – this is the 21st century, after all – is any media coverage and publicity, such as newspaper interviews, articles, and the like.

# Additional information
Running the following command to generante all explained project structure
```bash
  @prompt$: ./struct_project_generator.sh full_path_name_project 
```

# References
All exposed information is based on this page [http://nikola.me/folder_structure.html](http://nikola.me/folder_structure.html).