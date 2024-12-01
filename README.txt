1. Introduction
GRACE was designed to accept user input of target melting temperature (Tm) and specificity (deltaTm) and generate sequences of heterotrimeric collagen mimetic peptides. The generated heterotrimers are expected to form triple helices with melting temeperatures and specificities equal to or larger than user inputs. An additional feature of GRACE accepts user inputs of collagen-like peptide sequences and generate heterotrimer that features the inputs.

2. System requirements.
	a. macOS
	i.  GRACE was compiled using Xcode version 2408, thereby requiring macOS sonoma 14.5 or later to execute. 
 	b. <Linux version to be tested>

3. Installation guide 
	a. Download: Download the zipped GRACE file from SI [link].
	b.Unzip: Extract the contents of the ZIP file to a desired location. This will create a folder containing:
	i.GRACE: The executable file.
	ii. AminoAcids.csv: This file can be modified to limit the amino acid pool for sequence generation. Ensure that values are either 0 or 1.
	iii. parameters.txt: This file contains parameters for calculating Tm and deltaTm. Do not modify this file.


3. Executing GRACE
	a. Open the terminal 
	i. macOS: 
		Finder -> Applications -> Utilities -> Terminal 
	ii. Linux: <…>
	b. Navigate to the directory in which the unzipped GRACE is located 
		cd path/to/GRACE_directory
	c. Type command: 
		./GRACE
        d. Answer the prompts to start the search. 
	i. Searching time can range from minutes to hours depending on the given conditions and the length of sequences. 
	ii. Please expect that there are chances that the algorithm cannot find a set of sequences that satisfy all given conditions. In that case, the algorithm is terminated after 500,000 iterations and outputs the closet results found. 

4. Intepretation of outputs
	a. 



 
