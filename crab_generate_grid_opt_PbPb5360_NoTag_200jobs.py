import os

from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config


HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

config = config()

config.section_("General")
config.General.workArea = "crab_projects_grid"
config.General.requestName = "generate_grid_opt_PbPb5360_Notag_smoke_200jobs_v2"
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = "PrivateMC"
config.JobType.psetName = "upc_empty_source_cfg.py"
config.JobType.scriptExe = "run_generate_grid_opt_crab.sh"
config.JobType.inputFiles = [
    "run_generate_grid_opt_crab.sh",
    "upc_empty_source_cfg.py",
    "build_cmssw_1511/generate_grid_opt_job",
    "UPC_Probabilities_PbPb_5360GeV_Glauber_Final.root",
]
config.JobType.outputFiles = [
    "grid_PbPb5360_Notag_smoke_crab_job.tar",
]
config.JobType.allowUndistributedCMSSW = True
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 1500
config.JobType.maxJobRuntimeMin = 120

config.section_("Data")
# 400 M bins / 200 CRAB jobs = 2 M bins per job.  CRAB passes 1-based
# job ids to scriptExe; the wrapper converts them to generate_grid_opt_job's
# 0-based job_id.
config.Data.splitting = "EventBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = 200
config.Data.publication = False
config.Data.outputPrimaryDataset = "UPCGridOpt"
config.Data.outputDatasetTag = "generate_grid_opt_PbPb5360_Notag_smoke_200jobs_v2"
config.Data.outLFNDirBase = "/store/user/jianjie/RhoPrime4Pi/UPC/GridOpt_PbPb5360_Notag_smoke_200jobs_v2"

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"

# Jobs need only shipped input files.  If you want to force CERN, uncomment:
# config.Site.whitelist = ["T2_CH_CERN"]


if __name__ == "__main__":
    crabCommand("submit", config=config)
