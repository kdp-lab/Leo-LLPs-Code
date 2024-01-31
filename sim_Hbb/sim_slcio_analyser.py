import pyLCIO
import ROOT
import glob
import json
from math import *
import numpy as np


# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

# Gather input files
fnames = glob.glob("/local/d1/lrozanov/mucoll-tutorial-2023/sim_Hbb/stau_output_sim.slcio") 
print("Found %i files."%len(fnames))

# Create empty lists for each variable
mcp_pt = [] #mcp = MCParticle (truth)
mcp_phi = []
mcp_eta = []
pdgid = []
status = []

mcp_mu_pt = [] #mcp_mu = MCParticle muon
mcp_mu_phi = []
mcp_mu_eta = []

i = 0
# reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
# reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks", "SiTracks_Refitted", "MCParticle_SiTracks", "MCParticle_SiTracks_Refitted", "IBTrackerHits", "IETrackerHits", "OBTrackerHits", "OETrackerHits", "VBTrackerHits", "VETrackerHits"])
# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
for f in fnames:
    if f == "/data/fmeloni/DataMuC_MuColl10_v0A/reco/merged/muonGun_pT_1000_5000.slcio":
        continue
    if max_events > 0 and i >= max_events: break
    #if no_inner_hits > np.abs(max_events): break
    # if f != "/data/fmeloni/DataMuC_MuColl10_v0A/reco/merged/muonGun_pT_250_1000.slcio":
    #     print("Skipped")
    #     continue
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)
    for event in reader: 
        if max_events > 0 and i >= max_events: break
        #if no_inner_hits > np.abs(max_events): break
        if i%1 == 0: print("Processing event %i."%i)

        mcpCollection = event.getCollection("MCParticle")

        # Make counter variables
        n_mcp_mu = 0

        # MCPs
        imcp_pt = []
        imcp_eta = []
        imcp_phi = []
        ipdgid = []
        istatus = []

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_p = mcp.getMomentum()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
            pdg = mcp.getPDG()

            imcp_pt.append(mcp_tlv.Perp())
            imcp_eta.append(mcp_tlv.Eta())
            imcp_phi.append(mcp_tlv.Phi())
            ipdgid.append(pdg)
            istatus.append(mcp.getGeneratorStatus())
            # print("PID, Status, pT, eta, phi: ", pdg, mcp.getGeneratorStatus(), mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())

            if abs(mcp.getPDG())==13 and mcp.getGeneratorStatus()==1:

                imcp_mu_pt = []
                imcp_mu_eta = []
                imcp_mu_phi = []

                imcp_mu_pt.append(mcp_tlv.Perp())
                imcp_mu_eta.append(mcp_tlv.Eta())
                imcp_mu_phi.append(mcp_tlv.Phi())
                n_mcp_mu += 1
                # print("Truth pt, eta, phi:", mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())
                        
        # print("End of event \n")
        # This is here to check that we never reconstruct multiple muons
        # If we did, we'd have to match the correct muon to the MCP object to do eff/res plots
        # But since we don't, we can skip that step
        i += 1
        mcp_pt.append(imcp_pt)
        mcp_eta.append(imcp_eta)
        mcp_phi.append(imcp_phi)
        pdgid.append(ipdgid)
        status.append(istatus)
        if n_mcp_mu > 0:
            mcp_mu_pt.append(imcp_mu_pt)
            mcp_mu_eta.append(imcp_mu_eta)
            mcp_mu_phi.append(imcp_mu_phi)
       
    reader.close()

# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################
print("\nSummary statistics:")
print("Ran over %i events."%i)
print("Found:")
#print("\t%i MCPs"%hists["mcp_pt"].GetEntries())
#print("\t%i mu MCPs"%hists["mcp_mu_pt"].GetEntries())
print("\tSanity check mcpCollection length:", len(imcp_pt))

# Make a list of all the data you want to save
data_list = {}
data_list["mcp_pt"] = mcp_pt
data_list["mcp_eta"] = mcp_eta
data_list["mcp_phi"] = mcp_phi
data_list["pdgid"] = pdgid
data_list["status"] = status
data_list["mcp_mu_pt"] = mcp_mu_pt
data_list["mcp_mu_eta"] = mcp_mu_eta
data_list["mcp_mu_phi"] = mcp_mu_phi

# After the loop is finished, save the data_list to a .json file
output_json = "sim_Hbb/stau_output_sim.json"
with open(output_json, 'w') as fp:
    json.dump(data_list, fp)




