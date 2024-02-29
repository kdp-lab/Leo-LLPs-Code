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
fnames = glob.glob("/local/d1/lrozanov/mucoll-tutorial-2023/sim_Hbb/run_14_sim.slcio") 
print("Found %i files."%len(fnames))

# Create empty lists for each variable
mcp_pt = [] #mcp = monte-carlo particle (truth)
mcp_phi = []
mcp_eta = []
pdgid = []
status = []

VB_x, VB_y, VB_z, VB_time, VB_pdg = [], [], [], [], []
VE_x, VE_y, VE_z, VE_time, VE_pdg = [], [], [], [], []
IB_x, IB_y, IB_z, IB_time, IB_pdg = [], [], [], [], []
IE_x, IE_y, IE_z, IE_time, IE_pdg = [], [], [], [], []
OB_x, OB_y, OB_z, OB_time, OB_pdg = [], [], [], [], []
OE_x, OE_y, OE_z, OE_time, OE_pdg = [], [], [], [], []

i = 0
# reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
# reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks", "SiTracks_Refitted", "MCParticle_SiTracks", "MCParticle_SiTracks_Refitted", "IBTrackerHits", "IETrackerHits", "OBTrackerHits", "OETrackerHits", "VBTrackerHits", "VETrackerHits"])
# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
for f in fnames:
    if f == "/data/fmeloni/DataMuC_MuColl10_v0A/reco/merged/muonGun_pT_1000_5000.slcio":
        continue
    if max_events > 0 and i >= max_events: break
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)
    for event in reader: 
        if max_events > 0 and i >= max_events: break
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

         # Sim Hits
        iVB_x, iVB_y, iVB_z, iVB_time, iVB_pdg = [], [], [], [], []
        iVE_x, iVE_y, iVE_z, iVE_time, iVE_pdg = [], [], [], [], []
        iIB_x, iIB_y, iIB_z, iIB_time, iIB_pdg = [], [], [], [], []
        iIE_x, iIE_y, iIE_z, iIE_time, iIE_pdg = [], [], [], [], []
        iOB_x, iOB_y, iOB_z, iOB_time, iOB_pdg = [], [], [], [], []
        iOE_x, iOE_y, iOE_z, iOE_time, iOE_pdg = [], [], [], [], []

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

        for coll_name, (ix, iy, iz, itime, ipdg) in [
            ("VertexBarrelCollection", (iVB_x, iVB_y, iVB_z, iVB_time, iVB_pdg)),
            ("VertexEndcapCollection", (iVE_x, iVE_y, iVE_z, iVE_time, iVE_pdg)),
            ("InnerTrackerBarrelCollection", (iIB_x, iIB_y, iIB_z, iIB_time, iIB_pdg)),
            ("InnerTrackerEndcapCollection", (iIE_x, iIE_y, iIE_z, iIE_time, iIE_pdg)),
            ("OuterTrackerBarrelCollection", (iOB_x, iOB_y, iOB_z, iOB_time, iOB_pdg)),
            ("OuterTrackerEndcapCollection", (iOE_x, iOE_y, iOE_z, iOE_time, iOE_pdg))
        ]:
            try:
                collection = event.getCollection(coll_name)
                for hit in collection:
                    position = hit.getPosition()
                    ix.append(position[0])
                    iy.append(position[1])
                    iz.append(position[2])
                    itime.append(hit.getTime())

                    # Retrieve the MCParticle and its PDG code
                    mcp = hit.getMCParticle()
                    hit_pdg = mcp.getPDG() if mcp else None
                    ipdg.append(hit_pdg)

            except Exception as e:
                print(f"Error accessing {coll_name}: {e}")                
        # print("End of event \n")
        i += 1
        mcp_pt.append(imcp_pt)
        mcp_eta.append(imcp_eta)
        mcp_phi.append(imcp_phi)
        pdgid.append(ipdgid)
        status.append(istatus)
        VB_x.append(iVB_x); VB_y.append(iVB_y); VB_z.append(iVB_z); VB_time.append(iVB_time); VB_pdg.append(iVB_pdg)
        VE_x.append(iVE_x); VE_y.append(iVE_y); VE_z.append(iVE_z); VE_time.append(iVE_time); VE_pdg.append(iVE_pdg)
        IB_x.append(iIB_x); IB_y.append(iIB_y); IB_z.append(iIB_z); IB_time.append(iIB_time); IB_pdg.append(iIB_pdg)
        IE_x.append(iIE_x); IE_y.append(iIE_y); IE_z.append(iIE_z); IE_time.append(iIE_time); IE_pdg.append(iIE_pdg)
        OB_x.append(iOB_x); OB_y.append(iOB_y); OB_z.append(iOB_z); OB_time.append(iOB_time); OB_pdg.append(iOB_pdg)
        OE_x.append(iOE_x); OE_y.append(iOE_y); OE_z.append(iOE_z); OE_time.append(iOE_time); OE_pdg.append(iOE_pdg)


       
    reader.close()

# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################
print("\nSummary statistics:")
print("Ran over %i events."%i)
print("Found:")
print("\tSanity check mcpCollection length:", len(imcp_pt))

# Make a list of all the data you want to save
data_list = {
    "VB_x": VB_x, "VB_y": VB_y, "VB_z": VB_z, "VB_time": VB_time, "VB_pdg": VB_pdg,
    "VE_x": VE_x, "VE_y": VE_y, "VE_z": VE_z, "VE_time": VE_time, "VE_pdg": VE_pdg,
    "IB_x": IB_x, "IB_y": IB_y, "IB_z": IB_z, "IB_time": IB_time, "IB_pdg": IB_pdg,
    "IE_x": IE_x, "IE_y": IE_y, "IE_z": IE_z, "IE_time": IE_time, "IE_pdg": IE_pdg,
    "OB_x": OB_x, "OB_y": OB_y, "OB_z": OB_z, "OB_time": OB_time, "OB_pdg": OB_pdg
}
data_list["mcp_pt"] = mcp_pt
data_list["mcp_eta"] = mcp_eta
data_list["mcp_phi"] = mcp_phi
data_list["pdgid"] = pdgid
data_list["status"] = status

# After the loop is finished, save the data_list to a .json file
output_json = "run_14_sim.json"
with open(output_json, 'w') as fp:
    json.dump(data_list, fp)




