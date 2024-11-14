R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(libfastjet);
#include "fastjet/ClusterSequence.hh"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <iostream>
using namespace fastjet;
using namespace std;

void access_tree(){

//Event Class
erhic::EventMC *event(NULL); 

//Particle Class
erhic::ParticleMC *particle(NULL); 
erhic::ParticleMC *isitelectron(NULL); 
TFile *f = new TFile("../pythia8/outfiles/hepmcout.root");
TTree *tree = (TTree*)f->Get("EICTree");

Int_t nEntries = tree->GetEntries();
//Int_t nEntries =3;

cout<<"-------------------------------"<<endl;
cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer
//tree->SetBranchAddress("x", &x);
//tree->SetBranchAddress("Q2", &Q2);

Int_t nParticles;
Double_t Mult, Eta, Px, Py, Pz, Pt, E;
Int_t Status;
Int_t pdgCode;
Double_t charge;
Double_t totalJets;
int totalConstituents=0;
Int_t nPions = 0, nPionBars = 0, nKaons = 0, nKaonBars = 0;
Int_t nProtons = 0, nProtonBars = 0, nElectrons = 0, nPositrons = 0;
double x, Q2;
vector<PseudoJet> particles;
TLorentzVector part;
Int_t jetsize =0;
Int_t constituent_size;

//Single particle histograms
TH1D* hMult = new TH1D("hMult", "Multiplicity Distribution", 100, 0, 100);

TH1D* hEta1 = new TH1D("hEta1", "Charged Particle p_{T} Distribution", 100, 0, 10);
TH1D* hEta2 = new TH1D("hEta2", "Charged Particle p_{T} Distribution", 100, 0, 10);
TH1D* hEta3 = new TH1D("hEta3", "Charged Particle p_{T} Distribution", 100, 0, 10);

TH1D* hPion = new TH1D("hPion", "Pion Distribution", 100, 0, 1000);
TH1D* hPionBar = new TH1D("hPionBar", "Anti-Pion Distribution", 100, 0, 1000);
TH1D* hKaon = new TH1D("hKaon", "Kaon Distribution", 100, 0, 1000);
TH1D* hKaonBar = new TH1D("hKaonBar", "Anti-Kaon Distribution", 100, 0, 1000);
TH1D* hProton = new TH1D("hProton", "Proton Distribution", 100, 0, 1000);
TH1D* hProtonBar = new TH1D("hProtonBar", "Anti-Proton Distribution", 100, 0, 1000);
TH1D* hElectron = new TH1D("hElectron", "Electron Distribution", 100, 0, 1000);
TH1D* hPositron = new TH1D("hPositron", "Positron Distribution", 100, 0, 1000);

TH2D *hxVsQ2 = new TH2D("hxVsQ2", "2D Distribution of x Vs Q^{2}", 100, 0, 1, 100, 0, 100);

// histograms for jet properties
TH1D* hJetMultiplicity = new TH1D("hJetMultiplicity", "Jet Multiplicity;Number of Jet Constituents;Events", 10, 0, 10);
TH2D* hJetPtVsEta = new TH2D("hJetPtVsEta", "Jet p_{T} vs Jet #eta;Jet #eta;Jet p_{T} (GeV/c)", 100, -5, 5, 200, 0, 20);
TH2D* hJetPtVsQ2 = new TH2D("hJetPtVsQ2", "Q^{2} vs Jet p_{T};Jet p_{T} (GeV/c);Q^{2} ((GeV/c)^{2})", 100, 0, 100, 100, 0, 50);
TH2D* hJetPVsEta = new TH2D("hJetPVsEta", "Jet momentum vs Jet #eta;Jet #eta;Jet p (GeV/c)", 100, -5, 5, 200, 0, 20);
TH2D* hJetPVsQ2 = new TH2D("hJetPVsQ2", "Q^{2} vs Jet momentum;Jet p (GeV/c);Q^{2} ((GeV/c)^{2})", 100, 0, 100, 100, 0, 50);
TH1D* hJetsPerEvent = new TH1D("hNumJetsPerEvent", "Number of Jets per Event Distribution;Number of Jets;Events", 20, 0, 20);

//Loop Over Events
for(Int_t i=0;i<nEntries;i++)
{
    tree->GetEntry(i);
    Q2 = (Double_t) event->GetQ2();
    x = (Double_t) event->GetX();
    //printf("For Event %d, Q^2 = %.3f GeV^2!\n",i,Q2);
    hxVsQ2->Fill(x, Q2);
    nParticles = event->GetNTracks();
    printf("For Event %d, we have total %d particles!\n",i,nParticles);
    particles.clear();

    //Loop Over Each Particle
    for(Int_t j=0;j<nParticles;j++)
    {
        if (j==2) continue;
        particle = event->GetTrack(j);
        if (!particle) {
        std::cerr << "Warning: Null particle encountered at index " << j << ". Skipping...\n";
        continue;
        }
        //pdgCode = (Int_t) particle->GetPdgCode();
        Eta = (Double_t) particle->GetEta();
        Pt = (Double_t) particle->GetPt();
        Px = (Double_t) particle->GetPx();
        Py = (Double_t) particle->GetPy();
        Pz = (Double_t) particle->GetPz();
        E = (Double_t) particle->GetE();
        Status = (Int_t) particle->GetStatus(); 
        pdgCode = (Int_t) particle->Id();
        if (Status==1) Mult++;  
        if (Status != 1) continue;

        TDatabasePDG *pdgDatabase = TDatabasePDG::Instance(); 
        TParticlePDG *particlepdg = pdgDatabase->GetParticle(pdgCode); 

        if (particlepdg) {
            charge = particlepdg->Charge() / 3.0; // Charge is returned as 3*charge
            std::cout << "Charge of particle with PDG ID " << pdgCode << " is: " << charge << std::endl;
        } 

        if (abs(charge) == 1 || charge == 0)  // full jet selection
        //if (abs(charge) == 1)  // charged jet selection
        { 
            printf("For Event %d, particle %d: Id = %d\n",i,j,pdgCode);

            if (Status==1 && fabs(Eta)<3.5) particles.push_back(PseudoJet(Px, Py, Pz, E));  
            if (pdgCode == 211) {  // pi+
                nPions++;
            } else if (pdgCode == -211) {  // pi-
                nPionBars++;
            } else if (pdgCode == 321) {  // K+
                nKaons++;
            } else if (pdgCode == -321) {  // K-
                nKaonBars++;
            } else if (pdgCode == 2212) {  // p
                nProtons++;
            } else if (pdgCode == -2212) {  // anti-p
                nProtonBars++;
            } else if (pdgCode == 11) {  // e-
                nElectrons++;
            } else if (pdgCode == -11) {  // e+
                nPositrons++;
            }
        } // end charged/neutral selection loop
    }   //end particle loop

    // choose a jet definition, run the clustering, extract the jets
    double R = 1.0;
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    cout << "Clustering with " << jet_def.description() << endl;
    jetsize = jets.size();  

    for (unsigned ijet = 0; ijet < jets.size(); ijet++) 
    { 
        vector<PseudoJet> constituents = jets[ijet].constituents();
        constituent_size = constituents.size();
        totalConstituents += constituent_size;
        cout << "Event " << i << ", Jet " << ijet << ": Number of Constituents = " << constituent_size << endl;
        double jetPx = jets[ijet].px();
        double jetPy = jets[ijet].py();
        double jetPz = jets[ijet].pz();
        double jetP = sqrt(jetPx*jetPx + jetPy*jetPy + jetPz*jetPz);
        double jetPt = jets[ijet].pt();
        double jetEta = jets[ijet].eta();
        double jetPhi = jets[ijet].phi();    
        hJetPtVsEta->Fill(jetEta, jetPt);
        hJetPtVsQ2->Fill(Q2, jetPt);
        hJetPVsEta->Fill(jetEta, jetP);
        hJetPVsQ2->Fill(Q2, jetP);
        hJetMultiplicity->Fill(constituent_size);
        for (unsigned jjet = 0; jjet < constituents.size(); jjet++) {
            //hJetMultiplicity->Fill(constituent_size);
            //cout << "pT of constituents in each jet:" << constituents[j].pt() << endl;
            }
    }
    cout << "Number of jets in an event is" << jetsize <<  endl;
    hJetsPerEvent->Fill(jetsize); 
    totalJets += jetsize; 
}   //end event loop

Double_t meanJetsPerEvent = totalJets / nEntries;
Double_t meanConstituentsPerJet = totalConstituents / totalJets;

cout << "Number of pi+ (211): " << nPions << endl;
cout << "Number of pi- (-211): " << nPionBars << endl;
cout << "Number of K+ (321): " << nKaons << endl;
cout << "Number of K- (-321): " << nKaonBars << endl;
cout << "Number of protons (2212): " << nProtons << endl;
cout << "Number of anti-protons (-2212): " << nProtonBars << endl;
cout << "Number of electrons (11): " << nElectrons << endl;
cout << "Number of positrons (-11): " << nPositrons << endl;

cout << "Total number of particles produced is: " << Mult << endl;
cout << "Number of particles per event is: " << Mult/nEntries << endl;
cout << "Total number of jets in all events: " << totalJets << endl;
cout << "Mean number of jets per event: " << meanJetsPerEvent << endl;
cout << "Total number of constituents is: " << totalConstituents << endl;
cout << "Mean number of constituents per jet: " << meanConstituentsPerJet << endl;

hMult->Fill(Mult);
hPion->Fill(nPions);
hPionBar->Fill(nPionBars); 
hKaon->Fill(nKaons);
hKaonBar->Fill(nKaonBars); 
hProton->Fill(nProtons);
hProtonBar->Fill(nProtonBars); 
hElectron->Fill(nElectrons);
hPositron->Fill(nPositrons);

TCanvas *c1 = new TCanvas("c1", "2D Histogram of x vs Q^{2}", 800, 600);
c1->SetLogx(1);
c1->SetLogy(1);
c1->SetLogz(1);

hxVsQ2->Draw("COLZ");
hxVsQ2->SetStats(0);
hxVsQ2->GetXaxis()->SetTitle("x");
hxVsQ2->GetYaxis()->SetTitle("Q^{2} ((GeV/c)^{2})");

TPaveText *textBox = new TPaveText(0.35, 0.65, 0.65, 0.88, "NDC");
textBox->SetFillColor(kWhite);        // Set background color to white
textBox->SetBorderSize(1);            // Set border size
textBox->SetTextFont(42);             // Set text font
textBox->SetTextSize(0.04);           // Set text size
textBox->AddText(Form("Total Events: %d", nEntries));
textBox->AddText(Form("Total Particles: %f", Mult));
textBox->Draw();
c1->Update();

TCanvas *c2 = new TCanvas("c2", "Charged Particle p_{T}", 800, 600);
hEta1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
hEta1->GetYaxis()->SetTitle("Count");
hEta1->GetXaxis()->SetTitleOffset(0);
hEta1->GetYaxis()->SetTitleOffset(0);

hEta1->SetMarkerStyle(21);
hEta2->SetMarkerStyle(21);
hEta3->SetMarkerStyle(21);

hEta2->SetMarkerColor(kBlack);
hEta2->SetMarkerColor(kBlue);
hEta3->SetMarkerColor(kGreen);

hEta1->SetStats(0);
hEta2->SetStats(0);
hEta3->SetStats(0);

hEta1->Draw("p L E");
hEta2->Draw("p L E same");
hEta3->Draw("p L E same");
c2->Draw("L");

TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
leg->AddEntry(hEta1, "-3.5 < #eta < -1.5", "p");
leg->AddEntry(hEta2, "-1.5 < #eta < 1.5", "p");
leg->AddEntry(hEta3, "1.5 < #eta < 3.5", "p");
leg->Draw("same"); 

TLatex *jet = new TLatex();
jet->SetNDC();
jet->SetTextFont(42);
jet->SetTextSize(0.04);
jet->SetTextColor(kBlack);
jet->DrawLatex(0.3, 0.6, Form("Charged jet, anti-k_{T}, jet R =1.0"));
    
TFile* file = TFile::Open("output_fulljet_nov13.root", "RECREATE");
if (file && file->IsOpen()) {
    c1->Write("xVsQ2");
    c2->Write("pThistogram");
    hEta1->Write();
    hEta2->Write();
    hEta3->Write();
    hMult->Write();
    hPion->Write();
    hPionBar->Write();
    hKaon->Write();
    hKaonBar->Write();
    hProton->Write();
    hProtonBar->Write();
    hElectron->Write();
    hPositron->Write(); 
    hJetsPerEvent->Write();
    hJetMultiplicity->Write();
    hJetPtVsEta->Write();
    hJetPtVsQ2->Write();
    hJetPVsEta->Write();
    hJetPVsQ2->Write();
    file->Close();
}                                                                                         

}
