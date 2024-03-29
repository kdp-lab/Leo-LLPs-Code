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
max_events = -1 # Set to -1 to run over all events

# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
# fnames = glob.glob("/local/d1/lrozanov/mucoll-tutorial-2023/reco_Hbb_bib/134_0.1_reco_bib.slcio") 
fnames = glob.glob("/local/d1/mu+mu-/reco_bib/134_0.1_reco_bib.slcio") 
print("Found %i files."%len(fnames))

def check_hard_radiation(mcp, fractional_threshold):
    had_hard_rad = False
    daughters = mcp.getDaughters() 
    for d in daughters:
        if d.getPDG() == 22 or d.getPDG() == 23 or d.getPDG() == 24:        
            if d.getEnergy() > fractional_threshold*mcp.getEnergy():
                had_hard_rad = True   
    return had_hard_rad

# Create empty lists for each variable
mcp_pt = [] #mcp = MCParticle (truth)
mcp_phi = []
mcp_eta = []
pdgid = []
status = []
prod_vertex = []
prod_time = []
id = []

mcp_stau_pt = [] #mcp_stau = MCParticle stau
mcp_stau_phi = []
mcp_stau_eta = []

# TRACK
d0_res = [] #d0_res = d0 resolution
d0_res_match = [] #d0_res_match = d0 resolution for matched staus?
z0_res = [] #z0_res = z0 resolution
z0_res_match = [] #z0_res_match = z0 resolution for matched staus?

nhits = []
pixel_nhits = []
inner_nhits = []
outer_nhits = []
pt_res_hits = []

d0_res_vs_pt = [] 
d0_res_vs_eta = []
z0_res_vs_pt = []
z0_res_vs_eta = []
pt_res_vs_eta = [] #track muon pt resolution vs pt
pt_res_vs_pt = []
pt_res = [] 

# Truth matched 
pt_match = [] #This is truth pt
track_pt = [] #This is track pt
track_eta = [] #This is track eta
track_theta = []
eta_match = [] #This is truth eta
theta_match = []
phi_match = [] #This is track phi?
ndf = []
chi2 = []

# HITS
x = []
y = []
z = []
hit_pdgid = []
time = []
corrected_time = []
hit_layer = []
hit_detector = []
hit_side = []

sim_VB_x, sim_VB_y, sim_VB_z, sim_VB_time, sim_VB_pdg, sim_VB_mcpid, sim_VB_layer = [], [], [], [], [], [], []
sim_VE_x, sim_VE_y, sim_VE_z, sim_VE_time, sim_VE_pdg, sim_VE_mcpid, sim_VE_layer = [], [], [], [], [], [], []
sim_IB_x, sim_IB_y, sim_IB_z, sim_IB_time, sim_IB_pdg, sim_IB_mcpid, sim_IB_layer = [], [], [], [], [], [], []
sim_IE_x, sim_IE_y, sim_IE_z, sim_IE_time, sim_IE_pdg, sim_IE_mcpid, sim_IE_layer = [], [], [], [], [], [], []
sim_OB_x, sim_OB_y, sim_OB_z, sim_OB_time, sim_OB_pdg, sim_OB_mcpid, sim_OB_layer = [], [], [], [], [], [], []
sim_OE_x, sim_OE_y, sim_OE_z, sim_OE_time, sim_OE_pdg, sim_OE_mcpid, sim_OE_layer = [], [], [], [], [], [], []

reco_VB_x, reco_VB_y, reco_VB_z, reco_VB_time, reco_VB_layer = [], [], [], [], []
reco_VE_x, reco_VE_y, reco_VE_z, reco_VE_time, reco_VE_layer = [], [], [], [], []
reco_IB_x, reco_IB_y, reco_IB_z, reco_IB_time, reco_IB_layer = [], [], [], [], []
reco_IE_x, reco_IE_y, reco_IE_z, reco_IE_time, reco_IE_layer = [], [], [], [], []
reco_OB_x, reco_OB_y, reco_OB_z, reco_OB_time, reco_OB_layer = [], [], [], [], []
reco_OE_x, reco_OE_y, reco_OE_z, reco_OE_time, reco_OE_layer = [], [], [], [], []

speedoflight = 299792458/1000000  # mm/ns

# Make counter variables
n_mcp_stau = 0
no_inner_hits = 0
i = 0
num_matched_tracks = 0
num_dupes = 0
num_fake_tracks = 0
hard_rad_discard = 0
total_n_pfo_mu = 0
# reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
# reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks", "SiTracks_Refitted", "MCParticle_SiTracks", "MCParticle_SiTracks_Refitted", "IBTrackerHits", "IETrackerHits", "OBTrackerHits", "OETrackerHits", "VBTrackerHits", "VETrackerHits"])
# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
for f in fnames:
    if max_events > 0 and i >= max_events: break
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)
    for ievt,event in enumerate(reader): 
        if max_events > 0 and i >= max_events: break
        if i%10 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        trackCollection = event.getCollection("SiTracks")

        hit_collections = []
        IBTrackerHits = event.getCollection('ITBarrelHits')
        hit_collections.append(IBTrackerHits)
        IETrackerHits = event.getCollection('ITEndcapHits')
        hit_collections.append(IETrackerHits)
        OBTrackerHits = event.getCollection('OTBarrelHits')
        hit_collections.append(OBTrackerHits)
        OETrackerHits = event.getCollection('OTEndcapHits')
        hit_collections.append(OETrackerHits)
        VBTrackerHits = event.getCollection('VXDBarrelHits')
        hit_collections.append(VBTrackerHits)
        VETrackerHits = event.getCollection('VXDEndcapHits')
        hit_collections.append(VETrackerHits)
        
        # Relations
        VBrelationCollection = event.getCollection('VXDBarrelHitsRelations')
        VBrelation = pyLCIO.UTIL.LCRelationNavigator(VBrelationCollection)
        
        VErelationCollection = event.getCollection('VXDEndcapHitsRelations')
        VErelation = pyLCIO.UTIL.LCRelationNavigator(VErelationCollection)
        
        IBrelationCollection = event.getCollection('ITBarrelHitsRelations')
        IBrelation = pyLCIO.UTIL.LCRelationNavigator(IBrelationCollection)
        
        IErelationCollection = event.getCollection('ITEndcapHitsRelations')
        IErelation = pyLCIO.UTIL.LCRelationNavigator(IErelationCollection)
        
        OBrelationCollection = event.getCollection('OTBarrelHitsRelations')
        OBrelation = pyLCIO.UTIL.LCRelationNavigator(OBrelationCollection)
        
        OErelationCollection = event.getCollection('OTEndcapHitsRelations')
        OErelation = pyLCIO.UTIL.LCRelationNavigator(OErelationCollection)

        # MCPs
        imcp_pt = []
        imcp_eta = []
        imcp_phi = []
        ipdgid = []
        istatus = []
        iprod_vertex = []
        iprod_time = []
        iid = []

        # Tracks
        id0_res_vs_pt = []
        id0_res_vs_eta = []
        iz0_res_vs_pt = []
        iz0_res_vs_eta = []
        ipt_res_vs_eta = []
        ipt_res_vs_pt = []
        id0_res = []
        iz0_res = []
        ipt_res = []
        ipt_match = []
        itrack_pt = []
        itrack_eta = []
        itrack_theta = []
        ieta_match = []
        itheta_match = []
        iphi_match = []
        indf = []
        ichi2 = []
        id0_res_match = []
        iz0_res_match = []
        inhits = []

        # HITS
        ix = []
        iy = []
        iz = []
        ihit_pdg = []
        itime = []
        icorrected_time = []
        ihit_layer = []
        ihit_detector = []
        ihit_side = []
        
        # Sim Hits
        isim_VB_x, isim_VB_y, isim_VB_z, isim_VB_time, isim_VB_pdg, isim_VB_mcpid, isim_VB_layer = [], [], [], [], [], [], []
        isim_VE_x, isim_VE_y, isim_VE_z, isim_VE_time, isim_VE_pdg, isim_VE_mcpid, isim_VE_layer = [], [], [], [], [], [], []
        isim_IB_x, isim_IB_y, isim_IB_z, isim_IB_time, isim_IB_pdg, isim_IB_mcpid, isim_IB_layer = [], [], [], [], [], [], []
        isim_IE_x, isim_IE_y, isim_IE_z, isim_IE_time, isim_IE_pdg, isim_IE_mcpid, isim_IE_layer = [], [], [], [], [], [], []
        isim_OB_x, isim_OB_y, isim_OB_z, isim_OB_time, isim_OB_pdg, isim_OB_mcpid, isim_OB_layer = [], [], [], [], [], [], []
        isim_OE_x, isim_OE_y, isim_OE_z, isim_OE_time, isim_OE_pdg, isim_OE_mcpid, isim_OE_layer = [], [], [], [], [], [], []

        ireco_VB_x, ireco_VB_y, ireco_VB_z, ireco_VB_time, ireco_VB_layer = [], [], [], [], []
        ireco_VE_x, ireco_VE_y, ireco_VE_z, ireco_VE_time, ireco_VE_layer = [], [], [], [], []
        ireco_IB_x, ireco_IB_y, ireco_IB_z, ireco_IB_time, ireco_IB_layer = [], [], [], [], []
        ireco_IE_x, ireco_IE_y, ireco_IE_z, ireco_IE_time, ireco_IE_layer = [], [], [], [], []
        ireco_OB_x, ireco_OB_y, ireco_OB_z, ireco_OB_time, ireco_OB_layer = [], [], [], [], []
        ireco_OE_x, ireco_OE_y, ireco_OE_z, ireco_OE_time, ireco_OE_layer = [], [], [], [], []

        for coll_name, relation, (ix, iy, iz, itime, ipdg, imcpid, ilayer), (irel_x, irel_y, irel_z, irel_time, irel_layer) in [
            ("VertexBarrelCollection", VBrelation, (isim_VB_x, isim_VB_y, isim_VB_z, isim_VB_time, isim_VB_pdg, isim_VB_mcpid, isim_VB_layer), (ireco_VB_x, ireco_VB_y, ireco_VB_z, ireco_VB_time, ireco_VB_layer)),
            ("VertexEndcapCollection", VErelation, (isim_VE_x, isim_VE_y, isim_VE_z, isim_VE_time, isim_VE_pdg, isim_VE_mcpid, isim_VE_layer), (ireco_VE_x, ireco_VE_y, ireco_VE_z, ireco_VE_time, ireco_VE_layer)),
            ("InnerTrackerBarrelCollection", IBrelation, (isim_IB_x, isim_IB_y, isim_IB_z, isim_IB_time, isim_IB_pdg, isim_IB_mcpid, isim_IB_layer), (ireco_IB_x, ireco_IB_y, ireco_IB_z, ireco_IB_time, ireco_IB_layer)),
            ("InnerTrackerEndcapCollection", IErelation, (isim_IE_x, isim_IE_y, isim_IE_z, isim_IE_time, isim_IE_pdg, isim_IE_mcpid, isim_IE_layer), (ireco_IE_x, ireco_IE_y, ireco_IE_z, ireco_IE_time, ireco_IE_layer)),
            ("OuterTrackerBarrelCollection", OBrelation, (isim_OB_x, isim_OB_y, isim_OB_z, isim_OB_time, isim_OB_pdg, isim_OB_mcpid, isim_OB_layer), (ireco_OB_x, ireco_OB_y, ireco_OB_z, ireco_OB_time, ireco_OB_layer)),
            ("OuterTrackerEndcapCollection", OErelation, (isim_OE_x, isim_OE_y, isim_OE_z, isim_OE_time, isim_OE_pdg, isim_OE_mcpid, isim_OE_layer), (ireco_OE_x, ireco_OE_y, ireco_OE_z, ireco_OE_time, ireco_OE_layer))
        ]:
            try:
                collection = event.getCollection(coll_name)
                for hit in collection:
                    # print("hit",collection.getTypeName(), collection.getParameters())
                    position = hit.getPosition()
                    ix.append(position[0])
                    iy.append(position[1])
                    iz.append(position[2])
                    itime.append(hit.getTime())

                    # Retrieve the MCParticle and its PDG code
                    mcp = hit.getMCParticle()
                    hit_pdg = mcp.getPDG() if mcp else None
                    mcpid = mcp.id() if mcp else None
                    ipdg.append(hit_pdg)
                    imcpid.append(mcpid)
                    
                    # Get layer
                    encoding = collection.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                    decoder = pyLCIO.UTIL.BitField64(encoding)
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    detector = decoder["system"].value()
                    layer = decoder['layer'].value()
                    side = decoder["side"].value()
                    
                    ilayer.append(layer)
                    
                    recohit = relation.getRelatedFromObjects(hit)
                    if len(recohit) > 0:
                        recohit = recohit[0]
                        rel_position = recohit.getPosition()
                        irel_x.append(rel_position[0])
                        irel_y.append(rel_position[1])
                        irel_z.append(rel_position[2])
                        irel_time.append(recohit.getTime())
                        irel_layer.append(layer)
                        
            except Exception as e:
                print(f"Error accessing {coll_name}: {e}")   # Storing this to use for matching in the next loop

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
            iprod_vertex.append([mcp.getVertex()[i] for i in range(3)])
            iprod_time.append(mcp.getTime())
            iid.append(mcp.id())
            # print("PID, Status, pT, eta, phi: ", pdg, mcp.getGeneratorStatus(), mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())

            if abs(mcp.getPDG())==1000015 or abs(mcp.getPDG())==2000015:

                imcp_stau_pt = []
                imcp_stau_eta = []
                imcp_stau_phi = []

                imcp_stau_pt.append(mcp_tlv.Perp())
                imcp_stau_eta.append(mcp_tlv.Eta())
                imcp_stau_phi.append(mcp_tlv.Phi())
                n_mcp_stau += 1
                # print("Truth pt, eta, phi:", mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())

        # Loop over the track objects
        for track in trackCollection:
            Bfield = 5 #T, 3.57 for legacy
            theta = np.pi/2- np.arctan(track.getTanLambda())
            phi = track.getPhi()
            eta = -np.log(np.tan(theta/2))
            pt  = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
            track_tlv = ROOT.TLorentzVector()
            track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
            # dr = mcp_tlv.DeltaR(track_tlv) # I don't think this works
            nhitz = track.getTrackerHits().size()

            d0 = track.getD0()
            z0 = track.getZ0()

            id0_res.append(d0)
            iz0_res.append(z0)
            inhits.append(nhitz)
            itrack_pt.append(pt)
            itrack_eta.append(eta)
            itrack_theta.append(theta)
            iphi_match.append(track.getPhi())

            indf.append(track.getNdf())
            ichi2.append(track.getChi2())
            # print("Reco pt, eta, phi, nhits, dr:", pt, eta, phi, nhitz, dr)

            for j, particle_pt in enumerate(imcp_pt):
                particle_eta = imcp_eta[j]
                particle_phi = imcp_phi[j]
                mcp_tlv = ROOT.TLorentzVector()
                mcp_tlv.SetPtEtaPhiE(particle_pt, particle_eta, particle_phi, 0)
                dr = mcp_tlv.DeltaR(track_tlv)
                if dr < 0.005:
                    #print(particle_pt, pt)
                    ptres = (particle_pt - pt) / particle_pt
                    #print(j, dr)
                    num_matched_tracks += 1

                    id0_res_vs_pt.append([particle_pt, d0])
                    id0_res_vs_eta.append([particle_eta, d0])
                    iz0_res_vs_pt.append([particle_pt, z0])
                    iz0_res_vs_eta.append([particle_eta, z0])
                    ipt_res_vs_eta.append([particle_eta, ptres])
                    ipt_res_vs_pt.append([particle_pt, ptres])
                    ipt_match.append(particle_pt)
                    ieta_match.append(particle_eta)
                    id0_res_match.append(d0)
                    iz0_res_match.append(z0)
                    ipt_res.append(ptres)
                    #itheta_match.append(mcp_stau_theta[j])

            # Hit per track
            iix = []
            iiy = []
            iiz = []
            iihit_pdg = []
            iitime = []
            iicorrected_time = []
            iihit_layer = []
            iihit_detector = []
            iihit_side = []
            
            pixel_nhit = 0
            inner_nhit = 0
            outer_nhit = 0
            for hit in track.getTrackerHits():
                    try:
                        mcp = hit.getMCParticle()
                        hit_pdg = mcp.getPDG()
                    except:
                        hit_pdg = 0
                # now decode hits
                    encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                    decoder = pyLCIO.UTIL.BitField64(encoding)
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    detector = decoder["system"].value()
                    layer = decoder['layer'].value()
                    side = decoder["side"].value()
                    if detector == 1 or detector == 2:
                        pixel_nhit += 1
                    if detector == 3 or detector == 4:
                        inner_nhit += 1
                    if detector == 5 or detector == 6:
                        outer_nhit += 1
                    position = hit.getPosition()
                    pos_x = position[0]
                    pos_y = position[1]
                    pos_z = position[2]

                    d = sqrt(position[0]*position[0] + position[1]
                     * position[1] + position[2]*position[2])
                    tof = d/speedoflight

                    resolution = 0.03
                    if detector > 2:
                        resolution = 0.06

                    corrected_t = hit.getTime()*(1.+ROOT.TRandom3(ievt).Gaus(0., resolution)) - tof
                    
                    iix.append(pos_x)
                    iiy.append(pos_y)
                    iiz.append(pos_z)
                    iihit_pdg.append(hit_pdg)
                    iitime.append(hit.getTime())
                    iicorrected_time.append(corrected_t)
                    iihit_detector.append(detector)
                    iihit_layer.append(layer)
                    iihit_side.append(side)


            pixel_nhits.append([pixel_nhit])
            inner_nhits.append([inner_nhit])
            outer_nhits.append([outer_nhit])
            ix.append(iix)
            iy.append(iiy)
            iz.append(iiz)
            ihit_pdg.append(iihit_pdg)
            itime.append(iitime)
            icorrected_time.append(iicorrected_time)
            ihit_detector.append(iihit_detector)
            ihit_layer.append(iihit_layer)
            ihit_side.append(iihit_side)
                        
        # print("End of event \n")
        # This is here to check that we never reconstruct multiple muons
        # If we did, we'd have to match the correct muon to the MCP object to do eff/res plots
        # But since we don't, we can skip that step
        i+=1
        d0_res.append(id0_res)
        z0_res.append(iz0_res)
        nhits.append(inhits)
        track_pt.append(itrack_pt)
        track_eta.append(itrack_eta)
        track_theta.append(itrack_theta)
        phi_match.append(iphi_match)
        ndf.append(indf)
        chi2.append(ichi2)
        x.append(ix)
        y.append(iy)
        z.append(iz)
        hit_pdgid.append(ihit_pdg)
        time.append(itime)
        corrected_time.append(icorrected_time)
        hit_layer.append(ihit_layer)
        hit_detector.append(ihit_detector)
        hit_side.append(ihit_side)
        if len(id0_res_vs_pt) > 0:
            #pt_res_hits.append(ipt_res_hits)
            d0_res_vs_pt.append(id0_res_vs_pt)
            d0_res_vs_eta.append(id0_res_vs_eta)
            z0_res_vs_pt.append(iz0_res_vs_pt)
            z0_res_vs_eta.append(iz0_res_vs_eta)
            pt_res_vs_eta.append(ipt_res_vs_eta)
            pt_res_vs_pt.append(ipt_res_vs_pt)
            pt_res.append(ipt_res)
            pt_match.append(ipt_match)
            eta_match.append(ieta_match)
            theta_match.append(itheta_match)
            d0_res_match.append(id0_res_match)
            z0_res_match.append(iz0_res_match)
        mcp_pt.append(imcp_pt)
        mcp_eta.append(imcp_eta)
        mcp_phi.append(imcp_phi)
        pdgid.append(ipdgid)
        status.append(istatus)
        prod_vertex.append(iprod_vertex)
        prod_time.append(iprod_time)
        id.append(iid)
        if n_mcp_stau > 0:
            mcp_stau_pt.append(imcp_stau_pt)
            mcp_stau_eta.append(imcp_stau_eta)
            mcp_stau_phi.append(imcp_stau_phi)
        sim_VB_x.append(isim_VB_x); sim_VB_y.append(isim_VB_y); sim_VB_z.append(isim_VB_z); sim_VB_time.append(isim_VB_time); sim_VB_pdg.append(isim_VB_pdg); sim_VB_mcpid.append(isim_VB_mcpid); sim_VB_layer.append(isim_VB_layer)
        sim_VE_x.append(isim_VE_x); sim_VE_y.append(isim_VE_y); sim_VE_z.append(isim_VE_z); sim_VE_time.append(isim_VE_time); sim_VE_pdg.append(isim_VE_pdg); sim_VE_mcpid.append(isim_VE_mcpid); sim_VE_layer.append(isim_VE_layer)
        sim_IB_x.append(isim_IB_x); sim_IB_y.append(isim_IB_y); sim_IB_z.append(isim_IB_z); sim_IB_time.append(isim_IB_time); sim_IB_pdg.append(isim_IB_pdg); sim_IB_mcpid.append(isim_IB_mcpid); sim_IB_layer.append(isim_IB_layer)
        sim_IE_x.append(isim_IE_x); sim_IE_y.append(isim_IE_y); sim_IE_z.append(isim_IE_z); sim_IE_time.append(isim_IE_time); sim_IE_pdg.append(isim_IE_pdg); sim_IE_mcpid.append(isim_IE_mcpid); sim_IE_layer.append(isim_IE_layer)
        sim_OB_x.append(isim_OB_x); sim_OB_y.append(isim_OB_y); sim_OB_z.append(isim_OB_z); sim_OB_time.append(isim_OB_time); sim_OB_pdg.append(isim_OB_pdg); sim_OB_mcpid.append(isim_OB_mcpid); sim_OB_layer.append(isim_OB_layer)
        sim_OE_x.append(isim_OE_x); sim_OE_y.append(isim_OE_y); sim_OE_z.append(isim_OE_z); sim_OE_time.append(isim_OE_time); sim_OE_pdg.append(isim_OE_pdg); sim_OE_mcpid.append(isim_OE_mcpid); sim_OE_layer.append(isim_OE_layer)
        reco_VB_x.append(ireco_VB_x); reco_VB_y.append(ireco_VB_y); reco_VB_z.append(ireco_VB_z); reco_VB_time.append(ireco_VB_time); reco_VB_layer.append(ireco_VB_layer)
        reco_VE_x.append(ireco_VE_x); reco_VE_y.append(ireco_VE_y); reco_VE_z.append(ireco_VE_z); reco_VE_time.append(ireco_VE_time); reco_VE_layer.append(ireco_VE_layer)
        reco_IB_x.append(ireco_IB_x); reco_IB_y.append(ireco_IB_y); reco_IB_z.append(ireco_IB_z); reco_IB_time.append(ireco_IB_time); reco_IB_layer.append(ireco_IB_layer)
        reco_IE_x.append(ireco_IE_x); reco_IE_y.append(ireco_IE_y); reco_IE_z.append(ireco_IE_z); reco_IE_time.append(ireco_IE_time); reco_IE_layer.append(ireco_IE_layer)
        reco_OB_x.append(ireco_OB_x); reco_OB_y.append(ireco_OB_y); reco_OB_z.append(ireco_OB_z); reco_OB_time.append(ireco_OB_time); reco_OB_layer.append(ireco_OB_layer)
        reco_OE_x.append(ireco_OE_x); reco_OE_y.append(ireco_OE_y); reco_OE_z.append(ireco_OE_z); reco_OE_time.append(ireco_OE_time); reco_OE_layer.append(ireco_OE_layer)
       
    reader.close()

# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################
print("\nSummary statistics:")
print("Ran over %i events."%i)
print("Found:")
# print("\t%i MCPs"%len(np.ravel(mcp_pt)))
print("\t%i Stau MCPs"%n_mcp_stau)
print("\tSanity check mcpCollection length:", len(imcp_pt))
print("\tSanity check trackCollection length:", len(itrack_pt))
# print("\tSanity check # hits:", len(np.ravel(x)))
# print("\t%i PFOs"%hists["pfo_pt"].GetEntries())
# print("\t%i mu PFOs"%hists["pfo_mu_pt"].GetEntries())
# print('\t%i matched muon tracks'%(num_matched_tracks))
# print('\t%i duplicates eliminated'%num_dupes)
# print('\t%i hard radiations discarded'%hard_rad_discard)
# print('\t%i fake tracks'%num_fake_tracks)
# print('\t%i GeV'%np.max(mcp_stau_pt))


# Make a list of all the data you want to save
data_list = {
    "sim_VB_x": sim_VB_x, "sim_VB_y": sim_VB_y, "sim_VB_z": sim_VB_z, "sim_VB_time": sim_VB_time, "sim_VB_pdg": sim_VB_pdg, "sim_VB_mcpid": sim_VB_mcpid, "sim_VB_layer": sim_VB_layer,
    "sim_VE_x": sim_VE_x, "sim_VE_y": sim_VE_y, "sim_VE_z": sim_VE_z, "sim_VE_time": sim_VE_time, "sim_VE_pdg": sim_VE_pdg, "sim_VE_mcpid": sim_VE_mcpid, "sim_VE_layer": sim_VE_layer,
    "sim_IB_x": sim_IB_x, "sim_IB_y": sim_IB_y, "sim_IB_z": sim_IB_z, "sim_IB_time": sim_IB_time, "sim_IB_pdg": sim_IB_pdg, "sim_IB_mcpid": sim_IB_mcpid, "sim_IB_layer": sim_IB_layer,
    "sim_IE_x": sim_IE_x, "sim_IE_y": sim_IE_y, "sim_IE_z": sim_IE_z, "sim_IE_time": sim_IE_time, "sim_IE_pdg": sim_IE_pdg, "sim_IE_mcpid": sim_IE_mcpid, "sim_IE_layer": sim_IE_layer,
    "sim_OB_x": sim_OB_x, "sim_OB_y": sim_OB_y, "sim_OB_z": sim_OB_z, "sim_OB_time": sim_OB_time, "sim_OB_pdg": sim_OB_pdg, "sim_OB_mcpid": sim_OB_mcpid, "sim_OB_layer": sim_OB_layer,
    "sim_OE_x": sim_OE_x, "sim_OE_y": sim_OE_y, "sim_OE_z": sim_OE_z, "sim_OE_time": sim_OE_time, "sim_OE_pdg": sim_OE_pdg, "sim_OE_mcpid": sim_OE_mcpid, "sim_OE_layer": sim_OE_layer,
    "reco_VB_x": reco_VB_x, "reco_VB_y": reco_VB_y, "reco_VB_z": reco_VB_z, "reco_VB_time": reco_VB_time, "reco_VB_layer": reco_VB_layer,
    "reco_VE_x": reco_VE_x, "reco_VE_y": reco_VE_y, "reco_VE_z": reco_VE_z, "reco_VE_time": reco_VE_time, "reco_VE_layer": reco_VE_layer,
    "reco_IB_x": reco_IB_x, "reco_IB_y": reco_IB_y, "reco_IB_z": reco_IB_z, "reco_IB_time": reco_IB_time, "reco_IB_layer": reco_IB_layer,
    "reco_IE_x": reco_IE_x, "reco_IE_y": reco_IE_y, "reco_IE_z": reco_IE_z, "reco_IE_time": reco_IE_time, "reco_IE_layer": reco_IE_layer,
    "reco_OB_x": reco_OB_x, "reco_OB_y": reco_OB_y, "reco_OB_z": reco_OB_z, "reco_OB_time": reco_OB_time, "reco_OB_layer": reco_OB_layer,
    "reco_OE_x": reco_OE_x, "reco_OE_y": reco_OE_y, "reco_OE_z": reco_OE_z, "reco_OE_time": reco_OE_time, "reco_OE_layer": reco_OE_layer
}
data_list["mcp_pt"] = mcp_pt
data_list["mcp_eta"] = mcp_eta
data_list["mcp_phi"] = mcp_phi
data_list["pdgid"] = pdgid
data_list["status"] = status
data_list["prod_vertex"] = prod_vertex
data_list["prod_time"] = prod_time
data_list["id"] = id
data_list["mcp_stau_pt"] = mcp_stau_pt
data_list["mcp_stau_eta"] = mcp_stau_eta
data_list["mcp_stau_phi"] = mcp_stau_phi

data_list["d0_res"] = d0_res
data_list["z0_res"] = z0_res
data_list["nhits"] = nhits
data_list["pixel_nhits"] = pixel_nhits
data_list["inner_nhits"] = inner_nhits
data_list["outer_nhits"] = outer_nhits
data_list["pt_res_hits"] = pt_res_hits
data_list["d0_res_vs_pt"] = d0_res_vs_pt
data_list["d0_res_vs_eta"] = d0_res_vs_eta
data_list["z0_res_vs_pt"] = z0_res_vs_pt
data_list["z0_res_vs_eta"] = z0_res_vs_eta
data_list["pt_res_vs_eta"] = pt_res_vs_eta
data_list["pt_res_vs_pt"] = pt_res_vs_pt
data_list["pt_res"] = pt_res
data_list["pt_match"] = pt_match
data_list["track_pt"] = track_pt
data_list["track_eta"] = track_eta
data_list["track_theta"] = track_theta
data_list["eta_match"] = eta_match
data_list["theta_match"] = theta_match
data_list["phi_match"] = phi_match
data_list["ndf"] = ndf
data_list["chi2"] = chi2
data_list["d0_res_match"] = d0_res_match
data_list["z0_res_match"] = z0_res_match

data_list["x"] = x
data_list["y"] = y
data_list["z"] = z
data_list["hit_pdgid"] = hit_pdgid
data_list["time"] = time
data_list["corrected_time"] = corrected_time
data_list["hit_layer"] = hit_layer
data_list["hit_detector"] = hit_detector
data_list["hit_side"] = hit_side

# After the loop is finished, save the data_list to a .json file
output_json = "../reco_Hbb_bib/134_0.1_reco_bib.json"
with open(output_json, 'w') as fp:
    json.dump(data_list, fp)





