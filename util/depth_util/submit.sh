#!/bin/bash

sbatch \
	--job-name="depth" \
	--output="depth.%j.out" \
	--error="depth.%j.err" \
	driver.sbatch

