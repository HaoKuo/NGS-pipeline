

# set a destination directory for all job outputs
java -Dbackend.providers.Local.config.root=/set/new/destination/directory -jar cromwell.jar run ...

# set config file path for the Cromwell engine
java -Dconfig.file=/path/to/application.conf -jar cromwell.jar run ...

# set concurrent job numbers for each workflow
java -Dconcurrent-job-limit=5 -jar cromwell.jar run ...

# enable result caching in database
java -Dcall-caching.enabled=true -Dlookup-docker-hash=false -jar cromwell.jar run ...

# one can specify a options.json file to set the config terms of the Cromwell engine
java -jar cromwell.jar run wf.wdl wf_inputs.json wf_options.json

e.g.: options.json
{
    "final_workflow_outputs_dir" : "/Users/shlee/Desktop",
    "jes_gcs_root": "gs://my-bucket/workflows"
}

# Valid keys and their meanings:
## Global (use with any backend)
1. write_to_cache - Accepts values true or false. If false, the completed calls from this workflow will not be added to the cache. See the Call Caching section for more details.
2. read_from_cache - Accepts values true or false. If false, Cromwell will not search the cache when invoking a call (i.e. every call will be executed unconditionally). See the Call Caching section for more details.
3. final_workflow_log_dir - Specifies a path where per-workflow logs will be written. If this is not specified, per-workflow logs will not be copied out of the Cromwell workflow log temporary directory/path before they are deleted.
4. final_workflow_outputs_dir - Specifies a path where final workflow outputs will be written. If this is not specified, workflow outputs will not be copied out of the Cromwell workflow execution directory/path.
5. final_call_logs_dir - Specifies a path where final call logs will be written. If this is not specified, call logs will not be copied out of the Cromwell workflow execution directory/path.
6. default_runtime_attributes - A JSON object where the keys are runtime attributes and the values are defaults that will be used through the workflow invocation. Individual tasks can choose to override these values. See the runtime attributes section for more information.
7. continueOnReturnCode - Can accept a boolean value or a comma separated list of integers in a string. Defaults to false. If false, then only return code of 0 will be acceptable for a task invocation. If true, then any return code is valid. If the value is a list of comma-separated integers in a string, this is interpreted as the acceptable return codes for this task.
8. workflow_failure_mode - What happens after a task fails. Choose from:
  8.1 ContinueWhilePossible - continues to start and process calls in the workflow, as long as they did not depend on the failing call
  8.2 NoNewCalls - no new calls are started but existing calls are allowed to finish
  8.3 The default is NoNewCalls but this can be changed using the workflow-options.workflow-failure-mode configuration option.
9. backend - Override the default backend specified in the Cromwell configuration for this workflow only.
