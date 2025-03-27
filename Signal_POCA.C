#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <THStack.h>
#include <TMath.h>
#include <TStyle.h>

float dca_point_track(TVector3 pnt, TVector3 refv1, TVector3 p1) {
    // Calculate the DCA using the cross product and magnitude
    TVector3 p_tmp = pnt - refv1;
    float dca = p_tmp.Cross(p1).Mag() / p1.Mag();
    return dca;
}


int find_2ndvtx(TVector3 refv1, TVector3 p1, TVector3 refv2, TVector3 p2, TVector3 &vtx, float &dca) {
    TVector3 u(p1); // Momentum of particle 1
    TVector3 v(p2); // Momentum of particle 2
    TVector3 w = refv1 - refv2; // Vector connecting the two reference positions

    float a = u.Dot(u); // |u|^2
    float b = u.Dot(v); // u·v
    float c = v.Dot(v); // |v|^2
    float d = u.Dot(w); // u·w
    float e = v.Dot(w); // v·w
    float D = a * c - b * b; // Determinant for solving the system

    float sc, tc; // Parameters along momentum directions
    if (D < 1e-6) { // Handle degenerate case
        sc = 0;
        tc = (b > c ? d / b : e / c);
    } else {
        sc = (b * e - c * d) / D;
        tc = (a * e - b * d) / D;
    }

    // Calculate closest approach vector and midpoint
    TVector3 dP = w + (sc * u) - (tc * v); // Vector of closest approach
    TVector3 mP = 0.5 * (w + (sc * u) + (tc * v)); // Midpoint of closest approach

    // Secondary vertex and DCA
    vtx = refv2 + mP;
    dca = dP.Mag();
    //std::cout << "u.Mag(): " << u.Mag() << ", v.Mag(): " << v.Mag() << ", w.Mag(): " << w.Mag() << std::endl;
    //std::cout << "a: " << a << ", b: " << b << ", c: " << c << std::endl;
    //std::cout << "d: " << d << ", e: " << e << ", D: " << D << std::endl;
    return 0;
}

void Signal_POCA() {
   

    // Create TChains for LQ, NC, CC, and SIDIS events
   TChain *T_s = new TChain("events");
     // Load LQ SIGNAL datasets
  for (int i= 0; i <= 999000; i += 1000) {
      TString LQ5 = Form("/project/6041615/qunib/LQ/LQ5/recon_5/split_%d_to_%d_recon.root", i, i + 999);
      T_s->Add(LQ5);
  }


//for (int i = 1; i <= 5; ++i) {
  //  TString NC100 = Form(
    //    "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/24.12.0/epic_craterlake/DIS/NC/18x275/minQ2=100/pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_%d.*.eicrecon.tree.edm4eic.root",
      //  i);
   // T_s->Add(NC100);
   // }



     // Create histograms for each cut and process

    TH1F *h0_s = new TH1F("h0_s", "LQ Events vs Cuts", 13, 0, 13);    

     // Define cut thresholds
    const float cut_epzh = 18.0;
    const float cut_misspt_low = 1.0;
    const float cut_misspt_high = 9.0;
    const float cut_asso_deltaR = 1.0;
    const float cut_pion_pt = 1.0;
    const float cut_awayPt = 1.0; 
    const float cut_nearPt = 3.0;
    const float cut_3pion_pt = 3.0;
    const float cut_ave_dl = 0.05;
    const float cut_dR_sum = 0.3;
    const float cut_corMass = 2.0;
        
    // Create TTreeReaders for NC, CC, and SIDIS
    TTreeReader reader_s(T_s);


TTreeReaderArray<int> gFlavor_s(reader_s, "ReconstructedChargedParticles.PDG");
TTreeReaderArray<float> px_s(reader_s, "ReconstructedChargedParticles.momentum.x");
TTreeReaderArray<float> py_s(reader_s, "ReconstructedChargedParticles.momentum.y");
TTreeReaderArray<float> pz_s(reader_s, "ReconstructedChargedParticles.momentum.z");

TTreeReaderArray<float> vx_s(reader_s, "CentralTrackVertices.position.x");
TTreeReaderArray<float> vy_s(reader_s, "CentralTrackVertices.position.y");
TTreeReaderArray<float> vz_s(reader_s, "CentralTrackVertices.position.z");
TTreeReaderArray<int> prim_vtx_index_s(reader_s, "PrimaryVertices_objIdx.index");


TTreeReaderArray<float> rcTrkLoca2_s(reader_s, "CentralCKFTrackParameters.loc.a");
TTreeReaderArray<float> rcTrkLocb2_s(reader_s, "CentralCKFTrackParameters.loc.b");
TTreeReaderArray<float> rcTrkPhi2_s(reader_s, "CentralCKFTrackParameters.phi");



int maxEvents = 1000000; //Limit the number of events 
// ===========================
// Loop for LQ
// ===========================
int entryCount_s = 0;
while (reader_s.Next()) {
    if (entryCount_s >= maxEvents) break;
   //std::cout << "Event: " << entryCount_s << std::endl;
     entryCount_s++;

     h0_s->Fill(0);

  TVector3 vertex_rc(-999., -999., -999.);
  // Find Reconstructed primary vertex
  if (prim_vtx_index_s.GetSize() > 0) {  // Ensure there is a primary vertex
      int rc_vtx_index = prim_vtx_index_s[0];  // Get primary vertex index
      vertex_rc.SetXYZ(vx_s[rc_vtx_index], vy_s[rc_vtx_index], vz_s[rc_vtx_index]);
  }
  // Skip event if no valid vertex
   if (vertex_rc.X() == -999. || vertex_rc.Y() == -999. || vertex_rc.Z() == -999.) continue;

std::cout << "Reconstructed Primary Vertex: "
              << "X = " << vertex_rc.X() << " , "
              << "Y = " << vertex_rc.Y() << " , "
              << "Z = " << vertex_rc.Z() << " " << std::endl;

   // Fill event count for primary vertex selection
    vertex_rc.SetXYZ(vertex_rc.X(), vertex_rc.Y(), vertex_rc.Z());
    h0_s->Fill(1);


    TLorentzVector trk_lab, trk_headon; 
    TVector3 p3miss, p3miss2;
    TLorentzRotation lab2headon = TLorentzRotation().RotateY(25e-3).Boost(sin(25e-3), 0.0, 0.0);
    float spx = 0, spy = 0, spz = 0;
    float spx2 = 0, spy2 = 0, spz2 = 0;
    float Epzh2 = 0;

     // Loop over flavors and process each particle
    for (size_t j = 0; j < gFlavor_s.GetSize(); j++) {
	    int pdg =gFlavor_s[j];
        // Skip neutrinos and electrons
        if (fabs(pdg) == 12 || fabs(pdg) == 14 || fabs(pdg) == 16) continue;
        if (fabs(pdg) == 11) continue;  // Skip electrons

        // Assign mass based on particle flavor
        if (fabs(pdg) == 311) {
            trk_lab.SetXYZM(px_s[j], py_s[j], pz_s[j], 0.49368);  // Kaon mass
        } else if (fabs(pdg) == 2212) {
            trk_lab.SetXYZM(px_s[j], py_s[j], pz_s[j], 0.93827);  // Proton mass
        } else if (fabs(pdg) == 22) {
            trk_lab.SetXYZM(px_s[j], py_s[j], pz_s[j], 0);  // Photon mass
        } else {
            trk_lab.SetXYZM(px_s[j], py_s[j], pz_s[j], 0.13957);  // Pion mass
        }

        // Skip low momentum tracks
        if (trk_lab.Pt() < 0.001) continue;

        // Lab frame sum
        if (fabs(trk_lab.Eta()) < 3.0) {
            spx+= px_s[j];
            spy+= py_s[j];
            spz+= pz_s[j];

	}
         trk_headon = lab2headon * trk_lab;

            // Apply stricter kinematic cuts in the head-on frame
         if (fabs(trk_headon.Eta()) < 3.0 && trk_headon.Pt() > 0.1) {
             spx2+= trk_headon.Px();
             spy2+= trk_headon.Py();
             spz2+= trk_headon.Pz();
	     Epzh2+= (trk_headon.E() - trk_headon.Pz());
	 }
    }

    // Compute missing transverse momentum
    p3miss.SetXYZ(-1.*spx, -1.*spy, -1.*spz);
    p3miss2.SetXYZ(-1.*spx2, -1.*spy2, -1.*spz2);
    
    // Compute kinematic variables
    double pth2 = sqrt(spx2 * spx2 + spy2 * spy2);
    double y_jb = Epzh2 /2./18. ;  // Assuming 18 GeV beam energy
    double q2_jb = pth2 * pth2 / (1.0 - y_jb);

   // Apply selection cuts for Epzh and pt_miss
    if (Epzh2 < cut_epzh) continue;
    h0_s->Fill(2);

    if (p3miss.Pt() < cut_misspt_low || p3miss.Pt() > cut_misspt_high) continue;
    h0_s->Fill(3);

       
  
    //Loop and look for the 1st pion

    TLorentzVector trk_lab_1, trk_headon_1, trk_lab_2, trk_headon_2; 
    float nearPtSum = 0, awayPtSum = 0;
    float dcaV0_pi0;  
   
    for (size_t j = 0; j < gFlavor_s.GetSize(); j++) {
	int pdg_1 = gFlavor_s[j];
        if (pdg_1 != 211) continue;        
std::cout << "Selected PDG: " << pdg_1 << " at index " << j << std::endl;
    // Reconstructed closest approach (PCA)
    float d0_1 = rcTrkLoca2_s[j];  // loc.a
    float z0_1 = rcTrkLocb2_s[j];  // loc.b
    float az_1 = rcTrkPhi2_s[j];  // phi

    // Compute POCA
    float pcax_1 = -d0_1 * TMath::Sin(az_1);
    float pcay_1 =  d0_1 * TMath::Cos(az_1);
    float pcaz_1 =  z0_1;
std::cout << "POCA for pion " << j << ": "
              << "pcax = " << pcax_1 << ", "
              << "pcay = " << pcay_1 << ", "
              << "pcaz = " << pcaz_1 << std::endl;
   
     float dca2d_1 = TMath::Sqrt(TMath::Power(pcax_1 - vertex_rc.X(), 2) +
                               TMath::Power(pcay_1 - vertex_rc.Y(), 2));
     // Reco track momentum
    float px_1 = px_s[j];
    float py_1 = py_s[j];
    float pz_1 = pz_s[j];

    // Define reco track variables
    TVector3 ref_pnt(pcax_1, pcay_1, pcaz_1);
    TVector3 p3p(px_1, py_1, pz_1);

    TLorentzVector p4_pi;
    p4_pi.SetXYZM(px_1, py_1, pz_1, 0.13957);

    if (p3p.Pt() < cut_pion_pt) continue;

    float eta = p3p.Eta();
    float phi = p3p.Phi();

    trk_lab_1.SetXYZM(px_1, py_1, pz_1, 0.13957);
    trk_headon_1 = lab2headon * trk_lab_1;

    // Calculate DCA to primary vertex for the primary pion
     dcaV0_pi0 = dca_point_track(vertex_rc, ref_pnt, p3p);    

     int nPionAsso = 0;
     int nAwayTrk = 0;
     int nNearTrk = 0;
     TLorentzVector p4_3pi, p4_pi1, p4_pi2;
     TVector3  p3p_pion1, p3p_pion2;
     TVector3 ref_pnt_pion1, ref_pnt_pion2;
     float dca_pi01, dca_pi02, dca_pi12;
     float  dcaV0_pi01, dcaV0_pi02;
     

    //2nd loop to find near- and away- side tracks, find 2 associated pions
    for (size_t k = 0; k < gFlavor_s.GetSize(); k++) {
	    if (k == j) continue; // Skip the same track from the first loop

	 int pdg_2 = gFlavor_s[k];
         if (fabs(pdg_2) == 11) continue;  // Skip electrons

        if (pdg_2 != -211) continue;  // Only care about negative pions
std::cout << "Selected PDG: " << pdg_2 << " at index " << k << std::endl;
  	// Reconstructed closest approach (PCA)
	
    float d0_2 = rcTrkLoca2_s[k];  // loc.a
    float z0_2 = rcTrkLocb2_s[k];  // loc.b
    float az_2 = rcTrkPhi2_s[k];  // phi

    // Compute POCA
    float pcax_2 = -d0_2 * TMath::Sin(az_2);
    float pcay_2 =  d0_2 * TMath::Cos(az_2);
    float pcaz_2 =  z0_2;

        std::cout << "POCA for pion " << k << ": "
              << "pcax = " << pcax_2 << ", "
              << "pcay = " << pcay_2 << ", "
              << "pcaz = " << pcaz_2 << std::endl;

   float dca2d_2 = TMath::Sqrt(TMath::Power(pcax_2 - vertex_rc.X(), 2) +
                               TMath::Power(pcay_2 - vertex_rc.Y(), 2));
	
        // Reco track momentum
        float px_2 = px_s[k];
        float py_2 = py_s[k];
        float pz_2 = pz_s[k];

        TVector3 p3ptmp(px_2, py_2, pz_2);
	if (p3ptmp.Pt() < cut_pion_pt) continue;

         float eta_k = p3ptmp.Eta();
         float phi_k = p3ptmp.Phi();
                 
         trk_lab_2.SetXYZM(px_2, py_2, pz_2, 0.13957);
         trk_headon_2 = lab2headon*trk_lab_2;

         float deltaR = sqrt((eta - eta_k) * (eta - eta_k) + (phi - phi_k) * (phi - phi_k));
	 float deltaPhi = TVector2::Phi_mpi_pi(trk_headon_2.Phi() - trk_headon_1.Phi());
        //float DeltaPhi = trk_headon_2.DeltaPhi(trk_headon_1);
	 std::cout << "deltaR: " << deltaR << " deltaPhi: " << deltaPhi << std::endl;
	
         if (fabs(deltaPhi) < 2.0) {
            nearPtSum += trk_headon_2.Pt();
            nNearTrk++;
        }

        if ( fabs(trk_headon_2.DeltaPhi(-1.*trk_headon_1) ) < 2.0) {
            awayPtSum += trk_headon_2.Pt();
            nAwayTrk++;
	}

         //  if (pdg_2 != -211) continue;

            if (deltaR < cut_asso_deltaR) {
            nPionAsso++;
	    nNearTrk--;
	    nearPtSum -= trk_headon_2.Pt();

        }
       	else {
	 continue;
      	}
	
        if (nPionAsso == 1) {
           ref_pnt_pion1.SetXYZ(pcax_2, pcay_2, pcaz_2);
           p3p_pion1.SetXYZ(px_2, py_2, pz_2);
           p4_pi1.SetXYZM(px_2, py_2, pz_2, 0.13957);
           dcaV0_pi01 = dca_point_track(vertex_rc, ref_pnt_pion1, p3p_pion1);
        } else if (nPionAsso == 2) {
           ref_pnt_pion2.SetXYZ(pcax_2, pcay_2, pcaz_2);
           p3p_pion2.SetXYZ(px_2, py_2, pz_2);
           p4_pi2.SetXYZM(px_2, py_2, pz_2, 0.13957);
       	   dcaV0_pi02 = dca_point_track(vertex_rc, ref_pnt_pion2, p3p_pion2);
        }
} // done with second loop with k index


     // Only care about 3-pions 
        if (nPionAsso == 2) {
        h0_s->Fill(4);
        p4_3pi = p4_pi + p4_pi1 + p4_pi2;
  
       
        if (awayPtSum < cut_awayPt) continue;
        h0_s->Fill(5);

        if (nearPtSum > cut_nearPt) continue;
        h0_s->Fill(6);

        
        if (p4_3pi.Pt() < cut_3pion_pt) continue;
        h0_s->Fill(7);
	
         
        // Secondary vertex reconstruction
        TVector3 tau_vrt_candidate1, tau_vrt_candidate2, tau_vrt_candidate3;
        TVector3 tau_vrt_ave;

        find_2ndvtx(ref_pnt, p3p, ref_pnt_pion1, p3p_pion1, tau_vrt_candidate1, dca_pi01);
        find_2ndvtx(ref_pnt, p3p, ref_pnt_pion2, p3p_pion2, tau_vrt_candidate2, dca_pi02);
        find_2ndvtx(ref_pnt_pion1, p3p_pion1, ref_pnt_pion2, p3p_pion2, tau_vrt_candidate3, dca_pi12);
        

        float dl1 = (tau_vrt_candidate1 - vertex_rc).Mag();
        float dl2 = (tau_vrt_candidate2 - vertex_rc).Mag();
        float dl3 = (tau_vrt_candidate3 - vertex_rc).Mag();

        float dl_asy = fabs(dl1 - dl2) + fabs(dl1 - dl3) + fabs(dl2 - dl3);

        float ave_dl = (dl1 + dl2 + dl3) / 3.0;

        // Delta R calculations
        float dR_sum = (tau_vrt_candidate1 - vertex_rc).DeltaR(tau_vrt_candidate2 - vertex_rc) +
                       (tau_vrt_candidate1 - vertex_rc).DeltaR(tau_vrt_candidate3 - vertex_rc) +
                       (tau_vrt_candidate2 - vertex_rc).DeltaR(tau_vrt_candidate3 - vertex_rc);

        
        // Apply 30 µm vertex position requirement

        if (tau_vrt_candidate1.Mag() < 0.003 || tau_vrt_candidate2.Mag() < 0.003 || tau_vrt_candidate3.Mag() < 0.003) continue;
        h0_s->Fill(8);

	if(dl_asy > 0.05 ) continue;

        if (dR_sum > cut_dR_sum) continue;
        h0_s->Fill(9);

        if (ave_dl < cut_ave_dl) continue;
        h0_s->Fill(10);


        tau_vrt_ave = tau_vrt_candidate1 + tau_vrt_candidate2 + tau_vrt_candidate3;
        tau_vrt_ave = 1./3.*tau_vrt_ave;
        tau_vrt_ave = tau_vrt_ave - vertex_rc; 
                
        // Corrected mass cut
        float sintheta = p4_3pi.Vect().Cross(tau_vrt_ave).Mag() / tau_vrt_ave.Mag() / p4_3pi.Vect().Mag();
        float corMass = sqrt(p4_3pi.M2() + p4_3pi.P() * p4_3pi.P() * sintheta * sintheta) + p4_3pi.P() * sintheta;
        if (corMass > cut_corMass) continue;
        h0_s->Fill(11);

        // Final \(\Delta \phi\) cut
        if (p4_3pi.Vect().DeltaPhi(p3miss2) > 1.0) continue;
        h0_s->Fill(12);

      }
}//end of the loop where we started with 1st pion

}// end of while loop
   
THStack *hs = new THStack("hs", "");

hs->Add(h0_s);
h0_s->SetLineColor(kBlue+3);
h0_s->SetLineWidth(2);


TCanvas *c1 = new TCanvas("c1", "", 1400, 800);
 gPad->SetLogy(1);
 hs->Draw("HISTnostack");

hs->GetXaxis()->SetBinLabel(1, "input");
hs->GetXaxis()->SetBinLabel(2, "PrVtx");
hs->GetXaxis()->SetBinLabel(3, "Epzh");
hs->GetXaxis()->SetBinLabel(4, "misspt");
hs->GetXaxis()->SetBinLabel(5, "3#pi");
hs->GetXaxis()->SetBinLabel(6, "away1GeV");
hs->GetXaxis()->SetBinLabel(7, "nearIso");
hs->GetXaxis()->SetBinLabel(8, "3#pi_pt");
hs->GetXaxis()->SetBinLabel(9, "30#mum");
hs->GetXaxis()->SetBinLabel(10, "dRsum");
hs->GetXaxis()->SetBinLabel(11, "decayL");
hs->GetXaxis()->SetBinLabel(12, "cMass");
hs->GetXaxis()->SetBinLabel(13, "missing#Phi");

 TLegend *leg = new TLegend(0.61, 0.73, 0.81, 0.88);
  leg->AddEntry(h0_s, "#sigma_{V^{L}_{1/2}} = 4.7#times10^{-6}pb", "l");
  leg->SetTextSize(0.021);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->SaveAs("LQ5_Signal_RCP_POCA.png");
    
 }
