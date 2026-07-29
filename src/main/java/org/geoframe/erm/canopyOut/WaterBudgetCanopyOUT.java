/*
 * GNU GPL v3 License
 *
 * Copyright 2021 Niccolo` Tubini, Giuseppe Formetta
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
	
<<<<<<<< HEAD:src/main/java/it/geoframe/blogspot/canopy/WaterBudgetCanopyOUT.java
package it.geoframe.blogspot.canopy;

import static it.geoframe.blogspot.utility.Utils.getRKMean;
========
package org.geoframe.erm.canopyOut;

import static org.geoframe.erm.utility.Utils.getRKMean;
>>>>>>>> build/add-maven-pom:src/main/java/org/geoframe/erm/canopyOut/WaterBudgetCanopyOUT.java
import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;

import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;;

/**
 * The component solves the budget for the outer part of the canopy layer.
 * Inputs are: the rain and the potential evapotranspiration Outputs are: the
 * storage and the throughfall.
 * 
 * @author Marialaura Bancheri, Riccardo Busti, Daniele Andreis, Giuseppe Formetta
 * 
 */

public class WaterBudgetCanopyOUT {

	@Description("Input rain Hashmap")
	@In
	public HashMap<Integer, double[]> inHMRain;

	@Description("Input ETp Hashmap")
	@In
	public HashMap<Integer, double[]> inHMETp;

	@Description("Input CI Hashmap")
	@In
	public HashMap<Integer, double[]> initialConditionS_i;

	@Description("Leaf Area Index Hashmap")
	@In
	public HashMap<Integer, double[]> inHMLAI;

	@Description("coefficient canopy out")
	@In
	public double kc;

	@Description("Time step")
	@In
	public double tTimestep;

	@Description("Partitioning coefficient free throughfall")
	@In
	public double p;

	@Description("RK iterations")
	@In
	public double RKiter = 100;

	// @Description("ODE solver model:dp853, Eulero ")
	// @In
	// public String solver_model;

	@Description("The HashMap with the Actual input of the layer ")
	@Out
	public HashMap<Integer, double[]> outHMActualInput = new HashMap<Integer, double[]>();

	@Description("The HashMap with the Actual input of the layer ")
	@Out
	public HashMap<Integer, double[]> outHMActualOutput = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the Water Storage  ")
	@Out
	public HashMap<Integer, double[]> outHMStorage = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the Throughfall ")
	@Out
	public HashMap<Integer, double[]> outHMThroughfall = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the AET ")
	@Out
	public HashMap<Integer, double[]> outHMAET = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the AET ")
	@Out
	public HashMap<Integer, double[]> outHMError = new HashMap<Integer, double[]>();

	private HashMap<Integer, double[]>ciMap= new HashMap<Integer, double[]>();

	int step;
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
			double LAI = inHMLAI.get(ID)[0];
			double ETp = inHMETp.get(ID)[0];
			
			if (isNovalue(LAI))
				LAI = 0.6;
			LAI = (LAI == 0) ? 0.6 : LAI;
			if (step == 0) {
				double CI;
				if (initialConditionS_i != null) {
					CI = initialConditionS_i.get(ID)[0];
					if (isNovalue(CI))
						CI = kc * LAI / 2;
				} else {
					CI = kc * LAI / 2;
				}
				ciMap.put(ID, new double[] {CI});
			}
			double CI = ciMap.get(ID)[0];
			
			if (isNovalue(rain))
				rain = 0.0;
			if (isNovalue(ETp) || ETp < 0)
				ETp = 0.0;

			WaterBudgetCanopyStepResult r = calculateWaterBudgetCanopy(rain, LAI, ETp, CI, kc, p);
			
			// export to timeseries
			storeResult_series(ID, r.waterStorage, r.throughfall, r.AET, r.actualInput, r.actualOutput, r.error);

			// set new IC
			ciMap.put(ID, new double[] {r.waterStorage});
		}
		step++;
	}

	public static WaterBudgetCanopyStepResult calculateWaterBudgetCanopy(double rain, double LAI, double ETp, double CI, double kc, double p) {
		double s_CanopyMax = kc * LAI;

		double actualInput = (1 - p) * rain;

		// solve S at t^n+1
		double[] out = RK4(CI, actualInput, ETp, s_CanopyMax);
		double waterStorage = out[0];
		if (waterStorage < 0)
			waterStorage = 0;
		double error = out[1];

		// update variables at t^n+1
		double actualOutput = out[3];
		double throughfall = actualOutput + p * rain;
		double AET = out[2];
		
		WaterBudgetCanopyStepResult r = new WaterBudgetCanopyStepResult(
				waterStorage,
				throughfall,
				AET,
				actualInput,
				actualOutput,
				error
		);
		return r;
	}

	// compute dS/dt
	public static double[] computeFunction(double Sn, double in, double ETp, double s_CanopyMax) {
		if (Sn < 0) {
			Sn = 0;
		}
		double et = computeAET(Sn, in, ETp, s_CanopyMax);
		double actualOut = computeActualOutput(Sn, in, et, s_CanopyMax);
		return new double[] { in - et - actualOut, et, actualOut };
	}

	// compute AET
	public static double computeAET(double Sn, double in, double ETp, double s_CanopyMax) {
		return Math.min(Math.max(0, Sn + in), ETp * Math.min(1, Sn / s_CanopyMax));
	}

	// compute actual output
	public static double computeActualOutput(double Sn, double in, double et, double s_CanopyMax) {
		return Math.max(0, Sn + in - et - s_CanopyMax);
	}

	// RK4
	public static double[] RK4(double Sn, double in, double ETp, double s_CanopyMax) {

		double balance = 0;

		double[] k1 = computeFunction(Sn, in, ETp, s_CanopyMax);
		double[] k2 = computeFunction(Sn + 0.5 * k1[0], in, ETp, s_CanopyMax);
		double[] k3 = computeFunction(Sn + 0.5 * k2[0], in, ETp, s_CanopyMax);
		double[] k4 = computeFunction(Sn + k3[0], in, ETp, s_CanopyMax);
		double Sn1 = Sn + getRKMean(k1, k2, k3, k4, 0);
		double aet = getRKMean(k1, k2, k3, k4, 1);
		double actualOut = getRKMean(k1, k2, k3, k4, 2);
		;
		balance = balance + Sn - Sn1 + in - aet - actualOut;

		return new double[] { Sn1, balance, aet, actualOut };
	}

	// store results
	private void storeResult_series(int ID, double S, double tr, double aet, double in, double out, double err) {
		outHMStorage.put(ID, new double[] { S });
		outHMThroughfall.put(ID, new double[] { tr });
		outHMAET.put(ID, new double[] { aet });
		outHMActualInput.put(ID, new double[] { in });
		outHMActualOutput.put(ID, new double[] { out });
		outHMError.put(ID, new double[] { err });
	}
	
	public record WaterBudgetCanopyStepResult(//
			double waterStorage, //
			double throughfall, //
			double AET, //
			double actualInput, //
			double actualOutput, //
			double error //
	) {

	}

}