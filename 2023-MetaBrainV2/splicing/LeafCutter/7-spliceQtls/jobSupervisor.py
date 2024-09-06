import sys
import subprocess
import glob
import os
from pathlib import Path
from datetime import datetime
import time

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--dir", dest="jobdir",
	help="Directory to look for jobs", required=True)
parser.add_argument("--start", dest="jobstart",
	help="Job number to start with")
parser.add_argument("--stop", dest="jobstop",
	help="Job number to stop with")

args = vars(parser.parse_args())




result = subprocess.run(['whoami'], stdout=subprocess.PIPE)

user = result.stdout.decode('utf-8')
user = user.strip()
print("Running for user: "+user)

# cancel all running jobs
#print("Cancelling all running jobs")
#result = subprocess.run(['scancel','-u',user], stdout=subprocess.PIPE)

dt = datetime.now()
#print("{} - sleeping for 30 seconds while jobs are killed".format(dt))
#time.sleep(30)

maxjobs = 1000
sleeptime = 15
sleeptimeseconds = sleeptime * 60
devperc = 0
jobfolder = args["jobdir"]

jobstatus = {}

# to split jobs across clusters
# nb
jobStart = 0
jobStop = 1140

# gs
#jobStart = 1140
#jobStop = 2282

jobStart = -1
jobStop = -1

# get all jobs in folder
allJobs = []
for path in Path(jobfolder).rglob('*.sh'):
	allJobs.append(path)

allJobs.sort()
print("{} jobs found".format(len(allJobs)))

if jobStart > -1 and jobStop > -1:
	print("Focusing on jobs: {} - {}".format(jobStart, jobStop) )
	tmpJobs = []
	njobs = jobStop - jobStart
	for i in range(jobStart,jobStop):
		if i < len(allJobs):
			tmpJobs.append(allJobs[i])
	allJobs = tmpJobs

# determine which ones are finished
jobsRemain = []
jobsFinished = set()
for job in allJobs:
	finishedfile = str(job).replace('jobs','output')
	finishedfile = finishedfile.replace('.sh','-TopEffects.finished')
	if os.path.exists(finishedfile):
		jobsFinished.add(finishedfile)
	else:
		jobsRemain.append(job)

print("{} jobs finished, {} jobs remain".format(len(jobsFinished), len(jobsRemain)))

jobsRemain.sort()

# enter the while loop
while len(jobsRemain) > 0:
	# check the number of jobs in the queue
	result = subprocess.run(['squeue','-u',user], stdout=subprocess.PIPE)	
	lines = result.stdout.decode('utf-8').split("\n")
	queuectr = 0
	for line in lines:
		if len(line.strip()) > 0:
			queuectr += 1

	# if the number of jobs is smaller
	if queuectr < maxjobs:
		# submit some of the remaining jobs
		space = maxjobs - queuectr
		if space > len(jobsRemain):
			space = len(jobsRemain)

		# submit a set of jobs
		jctr = 0
		devctr = 0
		newJobRemain = []
		devjobs = round(space * devperc)
		print("{} - submitting {} jobs, {} dev jobs - {} remain in queue".format(dt, space, devjobs, len(jobsRemain)))
		for job in jobsRemain:
			if jctr < space + 1:
				if devctr < devjobs:
					result = subprocess.run(['sbatch','--qos=dev',job], stdout=subprocess.PIPE)
					devctr += 1
				else:
					result = subprocess.run(['sbatch',job], stdout=subprocess.PIPE)
				jctr += 1
			else:
				newJobRemain.append(job)
		jobsRemain = newJobRemain
		jobsRemain.sort()
		dt = datetime.now()
		print("{} - submitted {} jobs - {} remain in queue".format(dt, jctr, len(jobsRemain)))
	# sleep for 15 minutes
	dt = datetime.now()
	print("{} - sleeping for {} minutes ".format(dt, sleeptime))
	time.sleep(sleeptimeseconds)	
