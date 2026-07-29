/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

<<<<<<<< HEAD:src/main/java/it/geoframe/blogspot/rootzone/WaterBudgetRootZone.java
package it.geoframe.blogspot.rootzone;

import static it.geoframe.blogspot.utility.Utils.getRKMean;
========
package org.geoframe.erm.rootZone;

import static org.geoframe.erm.utility.Utils.getRKMean;
>>>>>>>> build/add-maven-pom:src/main/java/org/geoframe/erm/rootZone/WaterBudgetRootZone.java
import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

import java.util.HashMap;
import java.util.Set;
import java.util.Map.Entry;

import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;;
/**
 * The Class WaterBudget solves the water budget equation for the root zone
 * layer.
 * 
 * @author Marialaura Bancheri, Riccardo Busti, Daniele Andreis, Giuseppe Formetta
 */
public class WaterBudgetRootZone {

	@Description("Input rain Hashmap")
	@In
	public HashMap<Integer, double[]> inHMRain;

	@Description("Input ET wet canopy Hashmap")
	@In
	public HashMap<Integer, double[]> inHMEwc;

	@Description("Input ET Hashmap")
	@In
	public HashMap<Integer, double[]> inHMETp;

	@Description("Input CI Hashmap")
	@In
	public HashMap<Integer, double[]> initialConditionS_i;

	@Description("The maximum storage capacity")
	@In
	public double pCmax;

	@Description("Maximum percolation rate")
	@In
	public double g;

	@Description("Exponential of non-linear reservoir")
	@In
	public double h;

	@Description("Degree of spatial variability of the soil moisture capacity")
	@In
	public double pB_soil;

	// @Description("partitioning coefficient between the root zone and the runoff
	// reservoirs")
	// @Out
	// public double alpha;

	@Description("Maximum value of the water storage, needed for the computation of the Actual EvapoTraspiration")
	@In
	@Out
	public double s_RootZoneMax;

	@Description("CI of the water storage")
	@In
	@Out
	public double s_RootZoneCI;

	@Description("Initial saturation_degree")
	@In
	public double sat_degree = 0.5;

	@Description("RK iterations")
	@In
	public double RKiter = 100;

	// @Description("ODE solver model: dp853, Eulero")
	// @In
	// public String solver_model;

	@Description("The area of the HRUs in km2")
	@In
	public double A;

	@Description("Time step")
	@In
	public double tTimestep;

	@Description("The HashMap with the Actual input of the layer ")
	@Out
	public HashMap<Integer, double[]> outHMActualInput = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the Water Storage  ")
	@Out
	public HashMap<Integer, double[]> outHMStorage = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the AET ")
	@Out
	public HashMap<Integer, double[]> outHMEvaporation = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the outflow which drains to the lower layer")
	@Out
	public HashMap<Integer, double[]> outHMR = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the quick outflow ")
	@Out
	public HashMap<Integer, double[]> outHMquick = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the quick outflow ")
	@Out
	public HashMap<Integer, double[]> outHMquick_mm = new HashMap<Integer, double[]>();

	@Description("The output HashMap with alpha ")
	@Out
	public HashMap<Integer, double[]> outHMalpha = new HashMap<Integer, double[]>();

	@Description("The output HashMap with mass balance ")
	@Out
	public HashMap<Integer, double[]> outHMError = new HashMap<Integer, double[]>();

	int step;
	private HashMap<Integer, double[]>ciMap= new HashMap<Integer, double[]>();

	/**
	 * Process: reading of the data, computation of the storage and outflows
	 *
	 * @throws Exception the exception
	 */
	@Execute
	public void process() throws Exception {

		// reading the ID of all the stations
		Set<Entry<Integer, double[]>> entrySet = inHMRain.entrySet();

		// iterate over the station
		for (Entry<Integer, double[]> entry : entrySet) {
			Integer ID = entry.getKey();

			/** Input data reading */
			double rain = inHMRain.get(ID)[0];
			if (isNovalue(rain))
				rain = 0;

			double ETp = 0;
			if (inHMETp != null)
				ETp = inHMETp.get(ID)[0];
			if (isNovalue(ETp))
				ETp = 0;

			double Ewc = 0;
			if (inHMEwc != null)
				Ewc = inHMEwc.get(ID)[0];
			if (isNovalue(Ewc))
				Ewc = 0;

			if (step == 0) {
				double CI;
				if (initialConditionS_i != null) {
					CI = initialConditionS_i.get(ID)[0];
					if (isNovalue(CI))
						CI = s_RootZoneMax * sat_degree;

				} else {
					CI = s_RootZoneMax * sat_degree;
				}
				ciMap.put(ID, new double[] {CI});
			}
			double CI = ciMap.get(ID)[0];

			WaterBudgetRootZoneStepResult r = calculateWaterBudgetRootZone(rain, ETp, Ewc, CI, pB_soil, s_RootZoneMax, g,
					h, A, tTimestep);

			// export to timeseries
			storeResult_series(ID, r.waterStorage, r.actualInput, r.AET, r.recharge, r.quick_mm, r.quick, r.alpha, r.error);

			// set new IC
			ciMap.put(ID, new double[] {r.waterStorage});
		}
		step++;
	}

	public static  WaterBudgetRootZoneStepResult calculateWaterBudgetRootZone(double rain, double ETp, double Ewc, double CI, double pB_soil, double s_RootZoneMax, double g, double h, double A, double tTimestep) {
		double ETpNet = ETp - Ewc;
		double m3s = A * Math.pow(10, 3) / (tTimestep * 60);

		// solve S at t^n+1

		double[] out = RK4(CI, ETpNet, rain, pB_soil, s_RootZoneMax, g, h);
		double waterStorage = out[0];
		if (waterStorage < 0)
			waterStorage = 0;
		double error = out[1];

		// update variables at t^n+1
		double alfa = out[2];
		double quick_mm = out[3];
		double quick = quick_mm * m3s;
		double actualInput = out[4];
		double recharge = out[5];
		double AET = out[6];
		
		WaterBudgetRootZoneStepResult r = new WaterBudgetRootZoneStepResult(
				waterStorage,
				actualInput,
				AET,
				recharge,
				quick_mm,
				quick,
				alfa,
				error
		);

		// double alpha=(rain<0.001)?0:alpha(CI,rain,s_RootZoneMax);
		return r;
	}

	// compute dS/dt
	public static double[] computeFunction(double Sn, double etpnet, double rain, double pB_soil, double s_RootZoneMax, double g, double h) {
		if (Sn < 0) {
			Sn = 0;
		}
		double alpha = alpha(Sn, rain, pB_soil, s_RootZoneMax);
		double[] o = actualInputs(Sn, alpha, rain);
		double actualInputs = o[0];
		double quick = o[1];

		double aet = computeAET(Sn, actualInputs, etpnet, s_RootZoneMax);
		double recharge = computeR(Sn, actualInputs, aet, g, h, s_RootZoneMax);
		double fun = actualInputs - aet - recharge;
		return new double[] { fun, actualInputs, recharge, aet, alpha, quick };
	}

	// compute alpha according to Hymod
	private static double alpha(double Sn, double Pval, double pB_soil, double s_RootZoneMax) {
		double pCmax = s_RootZoneMax * (pB_soil + 1);
		double coeff1 = 1.0 - ((pB_soil + 1.0) * (Sn) / pCmax);
		double exp = 1.0 / (pB_soil + 1.0);
		double ct_prev = pCmax * (1.0 - Math.pow(coeff1, exp));
		double UT1 = Math.max((Pval - pCmax + ct_prev), 0.0);
		// Pval = Pval - UT1;
		double dummy = Math.min(((ct_prev + Pval - UT1) / pCmax), 1.0);
		double coeff2 = (1.0 - dummy);
		double exp2 = (pB_soil + 1.0);
		double xn = (pCmax / (pB_soil + 1.0)) * (1.0 - (Math.pow(coeff2, exp2)));
		double UT2 = Math.max(Pval - UT1 - (xn - Sn), 0);
		double alpha = (UT1 + UT2) / Pval;
		if (isNovalue(alpha) || alpha > 1)
			alpha = 1;
		// if (isNovalue(alpha)) alpha= 1;
		return alpha;
	}

	// compute actual inputs
	public static double[] actualInputs(double Sn, double alfa, double rain) {
		return new double[] { (1 - alfa) * rain, alfa * rain };
	}

	// compute groundwater recharge
	public static double computeR(double Sn, double in, double et, double g, double h, double s_RootZoneMax) {
		double out = g * Math.pow(Math.min(1, Sn / s_RootZoneMax), h);
		out = Math.min(Sn + in - et, out + Math.max(0, Sn - s_RootZoneMax + in - et - out));
		return out;

	}

	// compute AET
	public static double computeAET(double Sn, double in, double etpnet, double s_RootZoneMax) {
		return Math.min(Sn + in, etpnet * Math.min(1, 1.33 * Math.min(1, Sn / s_RootZoneMax)));
	}

	// RK4
	public static double[] RK4(double Sn, double etpnet, double rain, double pB_soil, double s_RootZoneMax, double g, double h) {

		double balance = 0;
		double[] k1 = computeFunction(Sn, etpnet, rain, pB_soil, s_RootZoneMax, g, h);
		double[] k2 = computeFunction(Sn + 0.5 * k1[0], etpnet, rain, pB_soil, s_RootZoneMax, g, h);
		double[] k3 = computeFunction(Sn + 0.5 * k2[0], etpnet, rain, pB_soil, s_RootZoneMax, g, h);
		double[] k4 = computeFunction(Sn + k3[0], etpnet, rain, pB_soil, s_RootZoneMax, g, h);
		double Sn1 = Sn + getRKMean(k1, k2, k3, k4, 0);
		double alpha = getRKMean(k1, k2, k3, k4, 4);
		double actualInput = getRKMean(k1, k2, k3, k4, 1);
		double quick = getRKMean(k1, k2, k3, k4, 5);
		double aet = getRKMean(k1, k2, k3, k4, 3);
		double recharge = getRKMean(k1, k2, k3, k4, 2);
		balance = balance + Sn - Sn1 + actualInput - aet - recharge;
		Sn = Sn1;
		return new double[] { Sn, balance, alpha, quick, actualInput, recharge, aet };
	}



	// store results
	private void storeResult_series(int ID, double S, double in, double aet, double re, double quick_mm, double quick,
			double alf, double err) {
		outHMStorage.put(ID, new double[] { S });
		outHMActualInput.put(ID, new double[] { in });
		outHMEvaporation.put(ID, new double[] { aet });
		outHMR.put(ID, new double[] { re });
		outHMquick_mm.put(ID, new double[] { quick_mm });
		outHMquick.put(ID, new double[] { quick });
		outHMalpha.put(ID, new double[] { alf });
		outHMError.put(ID, new double[] { err });

	}
	
	
	public record WaterBudgetRootZoneStepResult(
			double waterStorage,
			double actualInput,
			double AET,
			double recharge,
			double quick_mm,
			double quick,
			double alpha,
			double error
	) {
	}
	

}