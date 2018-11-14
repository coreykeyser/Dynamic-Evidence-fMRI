# Dynamic-Evidence-fMRI
Workflow basics:

If you haven't downloaded git on your computer, download it. And if you don't have a good code editor yet, I'd suggest trying atom. You should also download python if you haven't already....For the R coding, download R and Rtools and work through Rstudio.

Once you have that setup, go into your terminal and create a folder for project. Then cd into that folder it and type "git clone [the copied link from the green "clone" button in the top right of the main page]." This will download all of the materials to your computer. If you have any questions about how to run it on your computer tet corey at 3303280334.

https://education.github.com/git-cheat-sheet-education.pdf
HELPFUL GIT COMMANDS
to clone your directory-->"git clone [the copied link from the green "clone" button in the top right of the main page]"
to update your directory-->'git pull origin'
 
THE FILES IN THE GITHUB DO NOT INCLUDE THE EXPERIMENTAL SETUP SINCE THAT DOESN'T NEED MANIPULATED BUT IF YOU WANT TO BROWSE IT JUST REACH OUT TO ME AND I CAN GIVE YOU MORE DETAILED INSTRUCTIONS.
 
FILE GUIDE: the main data source is evidenceAccumulationMaster, this has the coherences for every trial along with the responses for 5 different subjects. "fast rlcca.r" runs a continuous versiin of the LCA. This is set up to take in a matrix of drift rates for each choice, here right or left. What is returned is a recording tracking of each accumulator through time. This is computed at each step in time by the parameter of the LCA, leak (K), lateral inhibition (L), and a noise term (eta). 

Main files:
DE Try 2.R
evidenceAccumulationMaster.csv	
fast rlcca.R
pseudolikelihood.r
realdata_track.R
runDe for DynEv.R
tracking.R
 
fMRI project:
creating a joint model between behvioral responses in LCA and the BOLD response from fMRI
 
Math Psych:
Strictly adapting, fitting, and performing parameter recovery on the continuous LCA....and maybe more stuff
1. adapting LCA code for the continuous experiment (done)
2. determining parameters - What does this mean? Fitting subjects?
3. performing parameter recovery
A. Workflow:
   Extract drift rates from a real experiment (I recommend only one trial)
   Run LCCA for many trials with the same set of drift rates and one set of parameters
   Estimate Leak and Lateral Inhibition
      Use you favourite MLE technique (its already set up for DE via Burnin)
         Use psuedoliklihood to approximate the likelihood of any potential parameters to fit the observed data (what you generated in step 2)
   Reapeat over other parameters of leak and lateral inhibition
   Try to reaoeat again but with non decision time
   Repeat with eta larger
   Repeat and try tonfind the least amount of data needed
   Sample from posterior using vanilla DEMCMC
         
   
