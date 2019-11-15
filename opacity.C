const Double_t cm2kpc = 3.2411945e-22;
TGraph *count_opacity(TGraph& dtaudx, Double_t distance) { // distance [kpc]
  TGraph *op = new TGraph();
  (*op) = dtaudx;
  for (Int_t i = 0; i < dtaudx.GetN(); i++)
    op->GetY()[i] = exp(- distance * dtaudx.GetY()[i] / cm2kpc);

  return op;
}

int opacity() {
  TGraph dtaudx("dtaudx_new.txt");

  const int NUM = 3;
  double Red[NUM]    = { 1.00, 0.00,  0.00 };
  double Green[NUM]  = { 0.00, 1.00,  0.00 };
  double Blue[NUM]   = { 0.00, 0.00,  1.00 };
  double Length[NUM] = { 0.00, 1.0/2, 1,   };

   const int nb = 10;
   int FI = TColor::CreateGradientColorTable(NUM,Length,Red,Green,Blue,nb);
   int MyPalette[nb];
   for (int i = 0; i < nb; i++) MyPalette[i]=FI+i;

  TCanvas can("opacity", "opacity", 800, 600);
  can.SetLogx();
  can.SetMargin(1.1, 0.03, 1.1, 0.03);
  TLegend leg(0.2, 0.2, 0.35, 0.6);
  TH1F *frame = can.DrawFrame(1e1, -0.05, 1e3, 1.05, ";Energy[TeV];e^{-#tau}");
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(0.81);
  frame->GetYaxis()->SetTitleSize(0.055);

  int iter = 0;
  for (Double_t distance = 100; distance <= 1000; distance += 100) {
    TGraph *op = count_opacity(dtaudx, distance);
    op->Draw("same l");
    op->SetLineColor(FI + iter);
    iter++;

    ostringstream title;
    title << distance << " kpc";
    leg.AddEntry(op, title.str().c_str(), "l");
  }

  TGaxis axis(1e1, 1.05, 1e3, 1.05, 1e1, 1e3, 510, "-G");
  axis.SetLabelSize(0);
  axis.Draw();

  TGaxis axisy(1e3, -0.05, 1e3, 1.05, -0.05, 1.05, 510, "+");
  axisy.SetLabelSize(0);
  axisy.Draw();

  leg.Draw();
  can.Print("opacity.eps");

  return 0;
}
