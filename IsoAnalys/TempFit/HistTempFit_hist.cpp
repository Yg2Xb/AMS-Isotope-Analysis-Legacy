#include "HistTempFit_hist.h"
#include "/afs/cern.ch/work/z/zuhao/public/yanzx/IsoAnalys/estimator/helper_func.cpp"
#include <TSystem.h>
#include <iostream>

using namespace AMS_Iso;

void HistTempFit(const std::string& isotype, const std::string& suffix,
		const std::string& inputDir, const std::string& outputDir, int UseMass) {

	// 获取同位素配置
	const auto& config = IsoFitConstants::configs.at(isotype);

	// 打开数据文件
	//TFile* f_data = TFile::Open((inputDir + "/" + config.dataName + "_temp_wide_bkg_L1Inner510cutoff.root ").c_str());
	//TFile* f_data = TFile::Open((inputDir + "/" + "B11_temp_wide_bkg_frag.root").c_str());
	TFile* f_data = TFile::Open((inputDir + "/" + "Bor_temp_narrow_bkg_unbL1_SDIATBeBcut.root").c_str());
	cout << f_data->GetName() << endl;
	//TFile* f_data2 = TFile::Open((inputDir + "/" + "B10_temp_wide_bkg_frag.root").c_str());

	// 打开MC文件
	std::vector<TFile*> f_mc;
	for (int mass : config.masses) {
		//std::string path = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/" + config.name + std::to_string(mass) + "_temp.root";
		//f_mc.push_back(TFile::Open(TString(path)));
		f_mc.push_back(TFile::Open((inputDir + "/" + config.name + std::to_string(mass) + "_temp_wide_bkgcut.root").c_str()));
	}

	// 创建画布
	TCanvas* canvas = new TCanvas("canvas", "Template Fit", 800, 600);
	TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
	canvas->Divide(1, 2);
	TPad* pad1 = (TPad*)canvas->cd(1);
	pad1->SetPad(0, 0.25, 1, 1);
	TPad* pad2 = (TPad*)canvas->cd(2);
	pad2->SetPad(0, 0, 1, 0.27);
	pad2->SetBottomMargin(0.3);
	pad2->SetGridy();

	// 创建PDF文件
	canvas->Print((outputDir + "/MassTF_" + config.name + "_" + suffix + 
				"wide_BtoBe_SDIATBeBcut_Use" + std::to_string(UseMass) + ".pdf[").c_str());

	// 获取能量区间
	const auto& ek_bins = Binning::WideBins;// getKineticEnergyBins(config.charge, UseMass);//
	for(int b = 0; b < ek_bins.size(); b++){
		cout<<ek_bins[b]<<endl;
	}
	// 创建用于存储所有探测器结果的ROOT文件
	TFile* output_file = TFile::Open(
			Form("%s/MassTF_%s_%swide_BtoBe_SDIATBeBcut_Use%d.root", outputDir.c_str(),
				config.name.c_str(), suffix.c_str(), UseMass),
			"RECREATE");

	const char* DetName[] = {"TOF", "NaF", "AGL"};
	int* beginBin[] = {};
	// 探测器循环
	for (int idet = 0; idet < 3; idet++) {
		// 设置拟合变量
		RooRealVar inv_mass("inv_mass", "1/mass", config.fitRangeLow[idet], 
				config.fitRangeUp[idet]);
		double alpha_min = 0.98, alpha_step = 0.002, alpha_max = 1.02;

		// 创建结果储存直方图
		TH2F* h_chi2 = new TH2F(Form("h_chi2_%s", DetName[idet]), 
				Form("%s Chi2 for each Ek bin", DetName[idet]), 
				ek_bins.size() - 1, ek_bins.data(), 21, alpha_min - alpha_step/2, 
				alpha_max + alpha_step/2);

		std::vector<TH2F*> h_fractions;
		std::vector<TH1F*> h_best_fractions;
		for (int i = 0; i < config.nFitParams; ++i) {
			h_fractions.push_back(new TH2F(
						Form("h_frac_%s%d_%s", config.name.c_str(), config.masses[i], 
							DetName[idet]),
						Form("%s %s%d fraction for each alpha and Ek bin", 
							DetName[idet], config.name.c_str(), config.masses[i]),
						ek_bins.size() - 1, ek_bins.data(), 21, alpha_min - alpha_step/2, 
						alpha_max + alpha_step/2));

			h_best_fractions.push_back(new TH1F(
						Form("h_best_%s%d_frac_%s", config.name.c_str(), 
							config.masses[i], DetName[idet]),
						Form("%s Best %s%d Fraction;Ek/n[GeV/n];%s%d Frac", 
							DetName[idet], config.name.c_str(), config.masses[i], 
							config.name.c_str(), config.masses[i]),
						ek_bins.size() - 1, ek_bins.data()));
		}

		// 创建其他结果直方图
		TH1F* h_best_alpha = new TH1F(Form("h_best_alpha_%s", DetName[idet]),
				Form("%s Best alpha;Ek/n[GeV/n];#alpha", DetName[idet]), ek_bins.size() - 1, ek_bins.data());
		TH1F* h_best_chi2 = new TH1F(Form("h_best_chi2_%s", DetName[idet]),
				Form("%s Best chi2/ndf;Ek/n[GeV/n];Chi2/NDF", DetName[idet]), ek_bins.size() - 1, ek_bins.data());
		TH1F* h_best_sgf = new TH1F(Form("h_best_sgf_%s", DetName[idet]),
				Form("%s Significance;Ek/n[GeV/n];Significance", DetName[idet]), ek_bins.size() - 1, ek_bins.data());
		TH1F* h_best_eff = new TH1F(Form("h_best_eff_%s", DetName[idet]),
				Form("%s Signal Eff;Ek/n[GeV/n];Signal Eff", DetName[idet]), ek_bins.size() - 1, ek_bins.data());
		TH1F* h_best_nsig = new TH1F(Form("h_best_nsig_%s", DetName[idet]),
				Form("%s Number of Signal;Ek/n[GeV/n];Number of Signal", DetName[idet]), ek_bins.size() - 1, ek_bins.data());
		TH1F* h_best_ntot = new TH1F(Form("h_best_ntot_%s", DetName[idet]),
				Form("%s Number of Total;Ek/n[GeV/n];Number of Total", DetName[idet]), ek_bins.size() - 1, ek_bins.data());
		TH1F* h_best_entries = new TH1F(Form("h_best_entries_%s", DetName[idet]),
				Form("%s Entries;Ek/n[GeV/n];Entries", DetName[idet]), ek_bins.size() - 1, ek_bins.data());

		// 获取数据直方图
		TH2F* data_hist_2d = (TH2F*)f_data->Get(Form("h_inv_%s_BorUseMass10",IsoFitConstants::x_titles[idet], config.dataName.c_str(), UseMass));
		data_hist_2d->RebinY(2);	
		//TH2F* data_hist_2d = (TH2F*)f_data->Get(Form("h_inv_%s_mass11_bin11_alpha10_FragToBe",IsoFitConstants::x_titles[idet]));
		//TH2F* data_hist_2d2 = (TH2F*)f_data2->Get(Form("h_inv_%s_mass10_bin10_alpha10_FragToBe",IsoFitConstants::x_titles[idet]));

		// 能量bin循环
		//for (int ibin = 6; ibin <= 40; ibin++) {
		for (int ibin = 1; ibin <= 14; ibin++) {
			// 检查bin是否在探测器范围内
			//if (idet == 0 && ibin >= config.detRanges[0]) continue;
			if (idet == 0 && ibin >= 4) continue;
			//if (idet == 1 && (ibin < 3 || ibin > 27)) continue;
			if (idet == 1 && (ibin <1 || ibin > 11)) continue;
			//if (idet == 2 && ibin < 9) continue;
			if (idet == 2 && ibin < 3) continue;
			//if(ibin > 30) continue;

			// 初始化本bin的最佳拟合记录
			FitCache fit_cache;
			double chi2_min = 5000;
			// alpha循环
			for (int j = 10; j < 11; j = j + 1) {
				double alpha = alpha_min + j * alpha_step;
				// 获取MC直方图
				std::vector<TH2F*> mc_hist_2d;
				for (int idm = 0; idm < config.masses.size(); idm++) {
					mc_hist_2d.push_back((TH2F*)f_mc[idm]->Get(Form("hL2original_inv_%s_mass%d_bin%d_alpha%d", 
									IsoFitConstants::x_titles[idet], config.masses[idm], UseMass, j)));
				}

				// 检查直方图有效性
				if (!data_hist_2d || std::any_of(mc_hist_2d.begin(), mc_hist_2d.end(),
							[](TH2F* h) { return !h; })) {
					continue;
				}

				// 投影数据直方图
				TH1D* data_hist = data_hist_2d->ProjectionX(Form("data_hist_bin%d_alpha%d", ibin, j), ibin+1, ibin+1);
				//TH1D* data_hist2 = data_hist_2d2->ProjectionX(Form("data2_hist_bin%d_alpha%d", ibin, j), ibin, ibin);
				//TH1D* data_hist = (TH1D*)data_hist1->Clone(Form("combined_hist_bin%d_alpha%d", ibin, j));
				//data_hist->Scale(0.7); 
				//data_hist->Add(data_hist2, 0.3); 
				data_hist->Rebin(2);
				data_hist->Sumw2();

				// 投影并平滑MC直方图
				std::vector<TH1D*> mc_hists;
				for (size_t i = 0; i < mc_hist_2d.size(); ++i) {
					TH1D* mc_hist = mc_hist_2d[i]->ProjectionX(
							Form("%s%d_UseMass%d_hist_bin%d_alpha%d", config.name.c_str(),
								config.masses[i], UseMass, ibin, j), ibin, ibin);
					cout<<"check temp binning, Bin10: "<<mc_hist_2d[i]->GetYaxis()->GetBinLowEdge(10)<<endl;
					cout<<"check temp Name: "<<mc_hist_2d[i]->GetName()<<endl;
					mc_hist->Rebin(2);
					mc_hist->Smooth(1);
					cout<<mc_hist->GetEntries()<<endl;;
					mc_hists.push_back(mc_hist);
				}

				// 检查直方图统计量
				//if (data_hist->GetMaximum() < 50 || data_hist->GetEntries() == 0 ||
				if (data_hist->GetMaximum() < 20 || data_hist->GetEntries() == 0 ||
						std::any_of(mc_hists.begin(), mc_hists.end(),
							[](TH1D* h) { return h->GetEntries() == 0; })) {
					// 清理本次循环的内存
					cout<<"data max:"<<data_hist->GetMaximum()<<endl;
					delete data_hist;
					for (auto* h : mc_hists) delete h;
					continue;
				}

				// 创建RooFit对象
				RooDataHist data("data", "Data", RooArgList(inv_mass), data_hist);
				std::vector<RooDataHist*> templates;
				std::vector<RooHistPdf*> pdfs;
				std::vector<RooRealVar*> fractions;

				// 创建模板和PDF
				for (size_t i = 0; i < mc_hists.size(); ++i) {
					templates.push_back(new RooDataHist(
								Form("%s%d_template", config.name.c_str(), config.masses[i]),
								Form("%s%d Template", config.name.c_str(), config.masses[i]),
								RooArgList(inv_mass), mc_hists[i]));

					pdfs.push_back(new RooHistPdf(
								Form("%s%d_pdf", config.name.c_str(), config.masses[i]),
								Form("%s%d PDF", config.name.c_str(), config.masses[i]),
								RooArgSet(inv_mass), *templates.back()));

					if (i < config.nFitParams) {
						fractions.push_back(new RooRealVar(
									Form("frac_%s%d", config.name.c_str(), config.masses[i]),
									Form("fraction of %s%d", config.name.c_str(), 
										config.masses[i]),
									config.initFractions[i], 0., 1.));
					}
				}

				// 创建剩余分数作为公式
				std::string formula = "1";
				RooArgList fracList;
				for (auto* frac : fractions) {
					formula += " - " + std::string(frac->GetName());
					fracList.add(*frac);
				}
				RooFormulaVar* last_frac = new RooFormulaVar(
						Form("frac_%s%d", config.name.c_str(), config.masses.back()),
						Form("fraction of %s%d", config.name.c_str(), config.masses.back()),
						formula.c_str(), fracList);

				// 创建模型
				RooArgList pdf_list;
				RooArgList frac_list;
				for (size_t i = 0; i < pdfs.size(); ++i) {
					pdf_list.add(*pdfs[i]);
					if (i < fractions.size()) {
						frac_list.add(*fractions[i]);
					} else {
						frac_list.add(*last_frac);
					}
				}
				RooAddPdf model("model", "Combined Model", pdf_list, frac_list);

				// 执行拟合
				cout<<"detector "<<DetName[idet]<<" begin fit: ek/n="<<ek_bins[ibin]<<" ,alpha= "<<alpha<<endl;
				RooFitResult* fit_res = model.fitTo(data, RooFit::Save(), RooFit::SumW2Error(kTRUE));
				cout<<"detector "<<DetName[idet]<<" begin fit: ek/n="<<ek_bins[ibin]<<" ,alpha= "<<alpha<<endl;


				RooPlot* frame = inv_mass.frame();
				data.plotOn(frame, RooFit::XErrorSize(0), RooFit::Name("data"));
				model.plotOn(frame, RooFit::LineColor(kRed), RooFit::LineWidth(3), 
						RooFit::Name("model"));

				double chi2 = calculateChi2(frame, "data", "model", config.fitRangeLow[idet], config.fitRangeUp[idet]);
				// 计算自由度
				int nBins = data_hist->FindBin(config.fitRangeUp[idet]) - data_hist->FindBin(config.fitRangeLow[idet]) + 1;
				int nParams = model.getParameters(data)->getSize();
				int ndf = nBins - nParams;
				double chi2_ndf = (ndf > 0 && fit_res->status() == 0) ? chi2/ndf : 5000;
				cout << "Chi2: " << chi2 << ", NDF: " << ndf << ", Chi2/NDF: " << chi2_ndf << endl;

				model.plotOn(frame, RooFit::VisualizeError(*fit_res, 1), 
						RooFit::FillColor(kRed), 
						RooFit::Name("model_with_errors"));

				const std::vector<int> lineColors = {kBlue, kMagenta, kOrange+1, kViolet};
				for (size_t i = 0; i < pdfs.size(); ++i) {
					model.plotOn(frame, RooFit::Components(*pdfs[i]), 
							RooFit::LineStyle(kSolid),
							RooFit::LineColor(lineColors[i % lineColors.size()]), 
							RooFit::LineWidth(3),
							RooFit::Name(Form("%s%d", config.name.c_str(), 
									config.masses[i])));
				}

				// 存储当前alpha的拟合结果到二维直方图
				h_chi2->SetBinContent(ibin, j + 1, chi2);
				for (size_t i = 0; i < fractions.size(); ++i) {
					h_fractions[i]->SetBinContent(ibin, j + 1, fractions[i]->getVal());
					h_fractions[i]->SetBinError(ibin, j + 1, fractions[i]->getError());
				}
				
				/*	
				SignificanceResult result = findOptimalSignificance(inv_mass, 
						pdfs.back(), data, model, *fractions[0], *fractions[1],
						config.fitRangeLow[idet], config.fitRangeUp[idet]);
				*/

				double signal_efficiency = 1;//result.efficiency;
				double upperBound = 1;//result.upperBound;
				double N_sig = 1;//result.N_sig;
				double N_tot = 1;//result.N_tot;
				double significance = 1;//result.significance;
				

				double n_entries = data_hist->Integral(data_hist->FindBin(config.fitRangeLow[idet]), data_hist->FindBin(config.fitRangeUp[idet]));
				//double n_entries = data_hist->Integral();

				// 如果是最佳拟合则更新缓存
				if (chi2_ndf < chi2_min) {

					fit_cache.UpdateCache(ibin, chi2_ndf, alpha, fit_res, fractions, frame, mc_hists, data_hist, significance, signal_efficiency, upperBound, N_sig, N_tot, n_entries);
					fit_cache.UpdateCache(ibin, chi2_ndf, alpha, fit_res, fractions, frame, mc_hists, data_hist, 1, 1, 1, 1, 1, n_entries);
					chi2_min = chi2_ndf;
				}

				// 清理本次循环的内存
				delete fit_res;
				delete frame;
				delete last_frac;
				for (auto* frac : fractions) delete frac;
				for (auto* pdf : pdfs) delete pdf;
				for (auto* temp : templates) delete temp;
				for (auto* h : mc_hists) delete h;
				delete data_hist;
			}  // alpha循环结束

			// 获取缓存中的最佳拟合结果
			const auto* best_fit = fit_cache.GetBestFit(ibin);
			if (best_fit && best_fit->is_valid) {

				// 绘制最佳拟合结果
				canvas->cd(1);
				best_fit->frame->SetXTitle(Form("1/%s", IsoFitConstants::x_titles[idet]));
				best_fit->frame->SetYTitle("Events");
				best_fit->frame->GetYaxis()->SetTitleOffset(1.2);
				best_fit->frame->SetTitle(Form("Ek/n in %.2fGeV/n - %.2fGeV/n", ek_bins[ibin - 1], ek_bins[ibin]));
				best_fit->frame->Draw();

				// 添加图例
				TLegend* legend = new TLegend(0.16, 0.61, 0.43, 0.85);
				legend->AddEntry(best_fit->frame->findObject("data"), "ISS Data", "ep");
				legend->AddEntry(best_fit->frame->findObject("model"), "Fitting", "l");
				for (size_t i = 0; i < config.masses.size(); ++i) {
					legend->AddEntry(best_fit->frame->findObject(
								Form("%s%d", config.name.c_str(), config.masses[i])),
							Form("%s%d Template", config.name.c_str(),
								config.masses[i]), "l");
				}
				setLegend(legend);
				legend->Draw("same");

				// 添加信息文本框
				TPaveText* pt = new TPaveText(0.71, 0.61, 0.85, 0.85, "NDC");
				pt->AddText(Form("Total ISS Entries: %.0f", best_fit->n_entries));
				for (size_t i = 0; i < best_fit->fractions.size(); ++i) {
					pt->AddText(Form("%s%d frac: %.3f #pm %.3f",
								config.name.c_str(), config.masses[i],
								best_fit->fractions[i],
								best_fit->fraction_errors[i]));
				}
				if (best_fit->is_valid) {
					pt->AddText(Form("Chi2/NDF: %.3f", best_fit->chi2_ndf));
					//pt->AddText(Form("Alpha: %.3f", best_fit->alpha));
				} else {
					pt->AddText("Fitting Failed");
				}
				/*
				pt->AddText(Form("%.2f Signal Efficiency:", best_fit->signal_efficiency));
				pt->AddText(Form("Significance(%.2f-%.2f): %.1f", config.fitRangeLow[idet], best_fit->upperBound, best_fit->significance));
				pt->AddText(Form("N_{sig}: %.1f", best_fit->n_signal));
				pt->AddText(Form("N_{tot}: %.1f", best_fit->n_total));
				*/
				setPaveText(pt);
				pt->Draw("same");

				// 绘制pull plot
				canvas->cd(2);
				TGraphErrors* pullGraph = new TGraphErrors();
				setupPullPlot(pullGraph, Form("1/%s",IsoFitConstants::x_titles[idet]), config.fitRangeLow[idet], config.fitRangeUp[idet]);
				calculatePull(best_fit->frame, pullGraph, "data", "model", config.fitRangeLow[idet], config.fitRangeUp[idet]);
				pullGraph->GetXaxis()->SetRangeUser(config.fitRangeLow[idet], config.fitRangeUp[idet]);
				pullGraph->GetYaxis()->SetRangeUser(-10, 10);
				pullGraph->Draw("AP");


				// 保存画布
				canvas->Print((outputDir + "/MassTF_" + config.name + "_" + suffix +
							"wide_BtoBe_SDIATBeBcut_Use" + std::to_string(UseMass) + ".pdf").c_str());

				// 更新结果直方图
				for (size_t i = 0; i < best_fit->fractions.size(); ++i) {
					h_best_fractions[i]->SetBinContent(ibin, best_fit->fractions[i]);
					h_best_fractions[i]->SetBinError(ibin, best_fit->fraction_errors[i]);
				}	
				h_best_alpha->SetBinContent(ibin, best_fit->alpha);
				h_best_chi2->SetBinContent(ibin, best_fit->chi2_ndf);
				h_best_sgf->SetBinContent(ibin, best_fit->significance);
				h_best_eff->SetBinContent(ibin, best_fit->signal_efficiency);
				h_best_nsig->SetBinContent(ibin, best_fit->n_signal);
				h_best_ntot->SetBinContent(ibin, best_fit->n_total);
				h_best_entries->SetBinContent(ibin, best_fit->n_entries);

				// 清理内存
				delete legend;
				delete pt;
				delete pullGraph;
			}
		}  // 能量bin循环结束

		// 将该探测器的所有结果写入文件
		output_file->cd();
		h_chi2->Write();
		for (auto* h : h_fractions) h->Write();
		for (auto* h : h_best_fractions) h->Write();
		h_best_alpha->Write();
		h_best_chi2->Write();
		h_best_sgf->Write();
		h_best_eff->Write();
		h_best_nsig->Write();
		h_best_ntot->Write();
		h_best_entries->Write();

		// 清理该探测器的内存
		delete h_chi2;
		for (auto* h : h_fractions) delete h;
		for (auto* h : h_best_fractions) delete h;
		delete h_best_alpha;
		delete h_best_chi2;
		delete h_best_sgf;
		delete h_best_eff;
		delete h_best_nsig;
		delete h_best_ntot;
		delete h_best_entries;
	}  // 探测器循环结束

	// 关闭文件
	output_file->Close();
	delete output_file;

	canvas->Print((outputDir + "/MassTF_" + config.name + "_" + suffix +
				"wide_BtoBe_SDIATBeBcut_Use" + std::to_string(UseMass) + ".pdf]").c_str());
	delete canvas;
	delete c2;

	f_data->Close();
	for (auto* f : f_mc) {
		f->Close();
		delete f;
	}
	delete f_data;
	}

	void HistTempFit_hist() {
		HistTempFit("Be", "", "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData",
				"/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit", 7);
		//HistTempFit("B", "", "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData",
		//		"/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit", 10);
		
		//   HistTempFit("Be", "", "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData",
		//   "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit", 9);
		//   HistTempFit("Be", "", "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData",
		//   "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit",10);
		   
	}
