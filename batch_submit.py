#!/usr/bin/env python
import argparse
import sys
import os
import datetime
from subprocess import Popen, PIPE, STDOUT
import logging

LOGGER_PREFIX = ' %s'
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def log(msg):
    logger.info(LOGGER_PREFIX % msg)




# -- CHANGE THESE

BOSON_MASS = 800

PT_HAT_MIN, PT_HAT_MAX = [200, 400]

# -- 

def simulation_dir():
    if os.environ['USER'] == 'lukedeo':
        return '/u/at/lukedeo/jet-simulations'
    elif os.environ['USER'] == 'bpn7':
        return '/nfs/slac/g/atlas/u01/users/bnachman/SLAC_pythia/Reclustering'
    else:
        raise ValueError('Invalid user!')


def generate_script(d):
    return 'cd {}\n'.format(simulation_dir()) + 'source ./setup.sh\n./event-gen/event-gen --OutFile {file} \
    --Proc {process} --NEvents {events} --pThatMin {pthatmin} --pThatMax {pthatmax} --BosonMass {bosonmass}'.format(**d)

def bsub_wrapper(script, name, queue, log):
    return 'bsub -o %s -q %s -J %s < %s' % (
            log, queue, name, script
        )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--jobs', type=int,
        help='Number of jobs to run', required=True)

    parser.add_argument('--events', type=int,
        help='Number of events', required=True)

    parser.add_argument('--output-dir', type=str, default='./files',
        help='path to directory where data files should be written')

    parser.add_argument('--file-prefix', type=str, default='GENERATED',
        help='prefix written to the start of all output files.')

    parser.add_argument('--log-dir', type=str, default='./logs',
        help='Path to directory to write the logs')
    
    parser.add_argument('--process', type=str, default='qcd',
        help='which process? One of {qcd, wprime}')


    process_dict = {'qcd' : 4, 'wprime' : 3}

    args = parser.parse_args()

    log('Determining output directories.')

    # -- get the simulation directory
    work_dir = simulation_dir()

    # -- create a timestamped subdir
    subdir = datetime.datetime.now().strftime('%b%d-%H%M%S')

    # -- make the logging directory
    if not os.path.exists(os.path.join(args.log_dir, subdir)):
        os.makedirs(os.path.join(args.log_dir, subdir))

    log('Will write logs to {}'.format(os.path.join(args.log_dir, subdir)))
    # -- make the job output directory
    scratch_space = os.path.join(args.output_dir, subdir)
    if not os.path.exists(os.path.join(args.output_dir, subdir)):
        os.makedirs(scratch_space)

    log('Will write samples to {}'.format(scratch_space))

    log('Script submitted for {} generation.'.format(args.process))
    # -- which process?
    process_code = process_dict[args.process]


    for job in xrange(args.jobs):
        log('Launching job %s of %s...' % (job + 1, args.jobs))
        log_file = os.path.join(args.log_dir, subdir, 'log_job_{}'.format(job))
        output_file = '%s_process_%s_bosonmass%s_pthat%s-%s_nevents%s_job%s.root' % (
                args.file_prefix, 
                args.process, 
                BOSON_MASS, 
                PT_HAT_MIN, 
                PT_HAT_MAX, 
                args.events,
                job
            )
        job_params = {
            'file' : os.path.join(scratch_space, output_file), 
            'process' : process_code, 
            'events' : args.events, 
            'pthatmin' : PT_HAT_MIN, 
            'pthatmax' : PT_HAT_MAX, 
            'bosonmass': BOSON_MASS
        }
        buf = os.path.join(work_dir, '.tmp-buf')
        with open(buf, 'wb') as _buf:
            _buf.write(generate_script(job_params))
        job_out = Popen(bsub_wrapper(buf, 'job_%s_of_%s' % (job, args.jobs), 'medium', log_file).split(), stdout=PIPE, stderr=STDOUT)
        log('Success.')

