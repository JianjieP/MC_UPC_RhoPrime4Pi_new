import FWCore.ParameterSet.Config as cms

process = cms.Process("UPCEMPTY")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))
