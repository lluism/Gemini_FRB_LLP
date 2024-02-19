## Gemini FRB Large and Long Program Code Repository

(Based on Charlie Kilpatrick's code in github https://github.com/charliekilpatrick/gemini_frb_llp)

Code, tools, and examples for planning and submitting observations, collecting data, and running pypeit in support of the Gemini FRB LLP (PI T. Eftekhari).


DOWNLOADING GEMINI DATA 

The important part is to obtain and save the cookie of the gemini search website. This cookie is used for downloading 
the data and identifying the user. 

to obtain the cookie (in Firefox and similar for Safari):
- open the gemini search website in the browser and log in. https://archive.gemini.edu/searchform
- in the website right-click and select Inspect, then Storage, and on the left the cookies will appear (archive.gemini.edu). 
- copy & paste the cookie and save it in a file; the cookie needs to be named like this: gemini_archive_session='my_archive_cookie' as detailed in the example of section 4.1 in https://archive.gemini.edu/help/api.html

Once the cookie is saved, follow the instructions in api/Example.md to run the .py file with your choice of data. 
Following the example this will save the downloaded data in a folder named /test.


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

RUNNING PYPEIT

The data saved in /test is split into folders but you should put everything in the same folder (e.g., /test) for pypeit 
to find it without issues later.

The script controlling the pypeit run is /pypeit/reduce_gmos.py

First enter the directory where all the observational data is (e.g., cd api/test), we will run pypeit from there. 



this is an example for running reduce_gmos.py:      python reduce_gmos.py ../api/test/ut220321/  gmos_north  

I have added -o string in the code to run_pypeit overwritting the previous runs






Known issues: 

2/10/2024: 

with the latest version of Pypeit, an error is raised because somewhere there is a np.float and that was 
deprecated in numpy v=1.20. I have downgraded numpy and corresponding packages to be able to run the code, but this 
needs to be corrected. 

Current Calibrations folder appears named Masters. I believe this is because I am using an old Pypeit version. 
The Science folder appears named science (lower case). I believe this is because I am using an old Pypeit version. 

