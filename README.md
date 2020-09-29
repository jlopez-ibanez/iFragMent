## 1. General Information

- This includes the program to calculate the score and p-value of a structural fingerprint (of a chemical entity) against the chemical profiles of its corresponding database following the methodology described in the manuscript by López-Ibáñez J., Pazos F. and Chagoyen M.

- The files necessary to run the calculations performed in the mentioned paper are included along with this program (excluding those from KEGG due to licencing problems).

- This program is written in python and makes use of ScIPy and NumPy packages so they must be properly installed.

## 2. Installation and use

- No installation is required as long as python and the libraries ScIPy and NumPy are available in your system.
	Information about how to install them can be found in the following URLs:
	
	- ScIPy and NumPy: https://scipy.org/install.html
	- Python: https://www.python.org/downloads/

- If neither python nor these packages are available in your system the most straightforward is to install them using a  distribution such as Anaconda, or its lightweight version miniconda (you may need to manually install Numpy and ScIPy):

	- Anaconda: https://www.anaconda.com/products/individual
	- Miniconda: https://docs.conda.io/en/latest/miniconda.html
 
- If executed without arguments, the program will run with default values. To change this behaviour it may be run via command line using as argument any of the IDs found in the ****entities**** folder (See *Section 3* for more information) as in the folowing example:
	> python ifragment.py c0044
- The program will calculate the score of the input compound against all the profiles of its corresponding database and save them in the ****results**** folder.

 
## 3. Files and Format

In the program's directory are included three folders:

##### ifragment_datafiles
Contains the files needed for the program to perform the calculations.  There are three different files for each database:
		
- Pathways file: the list of pairs pathways-compounds.
- Null model file: includes the pre-calculated parameters of the null model distribution for each of the pathways needed to calculate the z-score and p-value.
- Fingerprints file: a comma separated file including the compounds and structural descriptors of the database. Note that is not saved the whole fingerprints but only the index of the fragments.  

##### entities
Contains three files (one per database) with a list of the compound IDs and names for which can be calculated the scores. 

##### results
Files resulting from the execution of the program will be saved here preprending *scores_* to the compound ID passed as argument. Note that files will be overwritten without further advice.

## 4. Additional Notes

* You can calculate your own null model by randomizing a square matrix created from the fingerprints file.

* enviPath did not include any pathway ID so those included are not standard.
