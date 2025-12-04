
Setup:
* Download the data zip from the Boise State Scholarworks site.
* Unzip it so that there's a folder called 'data' at the same level as the 'code' folder. 'data/' should contain 'infrasound/', 'coordinates/', etc.
* Optionally change base_dir in the code files to the absolute path of the repo folder that contains 'code/', 'data/', etc. If you don't do this, everything will still work, but you'll just need to run all the code files from inside the 'code/' folder.

Running the code:
1. Run beamform.py to create the files in beamform_results/. Other code files require this output to exist before running them. Does not save any figures.
2. Run beam_stack.py to create the beam-stack spectrograms and save results. Does not save any figures.
3. Run the other files to analyze data and produce figures.
