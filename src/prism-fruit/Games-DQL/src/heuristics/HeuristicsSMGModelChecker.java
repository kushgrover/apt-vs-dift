//==============================================================================
//	
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package heuristics;

import java.util.HashSet;
import java.util.Set;

import heuristics.nextstate.HeuristicNextState;
import heuristics.nextstate.HeuristicNextStateFactory;
import heuristics.nextstate.HeuristicNextStateFactory.NextState;
import heuristics.search.Heuristic;
import heuristics.search.HeuristicFactory;
import heuristics.search.HeuristicFactory.HeuristicType;
import heuristics.search.HeuristicRTDP_Adj;
import heuristics.search.HeuristicSG;
import heuristics.update.StateUpdate;
import heuristics.update.StateUpdateFactory;
import heuristics.update.StateValueContainer;
import heuristics.update.ThreadSafeStateValueContainer;
import explicit.FastAdaptiveUniformisation;
import explicit.FastAdaptiveUniformisation.AnalysisType;
import explicit.ModelExplicit;
import explicit.ModelExplorer;
import explicit.ProbModelChecker.TermCrit;
import parser.Values;
import parser.ast.*;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismSettings;
import prism.Result;
import simulator.PrismModelExplorer;
import simulator.SimulatorEngine;

/**
 * SMG model checker based on heuristic methods.
 */
public class HeuristicsSMGModelChecker extends PrismComponent
{
	// Model file
	private ModulesFile modulesFile;
	// Properties file
	private PropertiesFile propertiesFile;
	// Simulator engine
	private SimulatorEngine engine;
	// Constants from model
	private Values constantValues;
	// Labels from the model
	private LabelList labelListModel;
	// Labels from the property file
	private LabelList labelListProp;
	//Heuristic name
	private HeuristicType type;
	//Next state type
	private NextState nextStateType;
	//Stopping criterion
	private double termCritParam;
	protected TermCrit termCrit = TermCrit.RELATIVE;
	private boolean verbose;
	private ExpressionTemporal currentExpression;
	private double epsilon;
	private int RTDP_ADJ_optimizations;
	private int statesOA;
	private int availableOA;
	private double delta;
	private double pmin;
	private int postMax;
	
	/**
	 * Constructor.
	 */
	public HeuristicsSMGModelChecker(PrismComponent parent, ModulesFile modulesFile, PropertiesFile propertiesFile, SimulatorEngine engine, String heuristicName, String nextStateName, boolean verbose, double epsilon) throws PrismException
	{
		super(parent);
		this.modulesFile = modulesFile;
		this.propertiesFile = propertiesFile;
		this.engine = engine;

		// Get combined constant values from model/properties
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		if (propertiesFile != null)
			constantValues.addValues(propertiesFile.getConstantValues());
		this.labelListModel = modulesFile.getLabelList();
		this.labelListProp = propertiesFile.getLabelList();
		this.type = HeuristicType.valueOf(heuristicName);
		this.nextStateType = NextState.valueOf(nextStateName);
		this.termCritParam = settings.getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM);
		this.verbose = verbose;
		this.epsilon = epsilon;
	}

	public HeuristicsSMGModelChecker(PrismComponent parent, ModulesFile modulesFile, PropertiesFile propertiesFile, SimulatorEngine engine, String heuristicName, String nextStateName, boolean verbose, double epsilon, int RTDP_ADJ_optimizations, int statesOA, int availableOA, double delta, double pmin, int post) throws PrismException
	{
		super(parent);
		this.modulesFile = modulesFile;
		this.propertiesFile = propertiesFile;
		this.engine = engine;

		// Get combined constant values from model/properties
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		if (propertiesFile != null)
			constantValues.addValues(propertiesFile.getConstantValues());
		this.labelListModel = modulesFile.getLabelList();
		this.labelListProp = propertiesFile.getLabelList();
		this.type = HeuristicType.valueOf(heuristicName);
		this.nextStateType = NextState.valueOf(nextStateName);
		this.termCritParam = settings.getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM);
		this.verbose = verbose;
		this.epsilon = epsilon;
		this.RTDP_ADJ_optimizations = RTDP_ADJ_optimizations;
		this.statesOA=statesOA;
		this.availableOA=availableOA;
		this.delta=delta;
		this.pmin=pmin;
		this.postMax=post;
	}
	
	public ModulesFile getModulesFile() {
		return modulesFile;
	}
	
	public PropertiesFile getPropertiesFile() {
		return propertiesFile;
	}
	
	public SimulatorEngine getEngine() {
		return engine;
	}
	
	public ExpressionTemporal getExpression() {
		return currentExpression;
	}

	/**
	 * Model check a property.
	 */
	public Result check(Expression expr) throws PrismException
	{
		Result res;
		String resultString;
		long timer;

		// Starting model checking
		timer = System.currentTimeMillis();

		// Do model checking
		res = checkExpression(expr);

		// Model checking complete
		timer = System.currentTimeMillis() - timer;
		mainLog.println("\nModel checking completed in " + (timer / 1000.0) + " secs.");

		// Print result to log
		resultString = "Result";
		if (!("Result".equals(expr.getResultName())))
			resultString += " (" + expr.getResultName().toLowerCase() + ")";
		resultString += ": " + res;
		mainLog.print("\n" + resultString + "\n");

		// Return result
		return res;
	}
	
	protected double getTermCritParam() {
		return termCritParam;
	}
	
	protected TermCrit getTermCrit() {
		return termCrit;
	}

	/**
	 * Model check an expression (used recursively).
	 */
	
	private ExpressionPATL lastExprRPATL;
	
	private Result checkExpression(Expression expr) throws PrismException
	{

		// Current range of supported properties is quite limited...
		if(expr instanceof ExpressionPATL) {
			ExpressionPATL exprRPATL = (ExpressionPATL)expr;
			lastExprRPATL = exprRPATL;
			if(exprRPATL.getExpressionType() == ExpressionPATL.PRB) {
				return checkExpressionProb((ExpressionProb) exprRPATL.getExpressionProb());
			}
			//else if (expr instanceof ExpressionReward)
				//res = checkExpressionReward((ExpressionReward) expr);
		}
		if(expr instanceof ExpressionProb) {
			ExpressionProb exprProb = (ExpressionProb)expr;
			return checkExpressionProb((ExpressionProb) exprProb);
			//else if (expr instanceof ExpressionReward)
				//res = checkExpressionReward((ExpressionReward) expr);
		}
		throw new PrismException("Only P queries are supported by heuristics");
	}

	/**
	 * Model check a P operator.
	 */
	private Result checkExpressionProb(ExpressionProb expr) throws PrismException
	{
		boolean min = false;
		// Check whether P=? (only case allowed)
		if (expr.getProb() != null) {
			throw new PrismException("Only P=? properties are supported");
		}

		if (!(expr.getExpression() instanceof ExpressionTemporal)) {
			throw new PrismException("Heuristic model checking currently only supports simple path operators");
		}
		ExpressionProb exprProb = (ExpressionProb)expr;
		if(exprProb.getRelOp().isMin()) {
			min = true;
		}
		
		ExpressionTemporal exprTemp = (ExpressionTemporal) expr.getExpression();
		currentExpression = exprTemp;
		if (!exprTemp.isSimplePathFormula()) {
			throw new PrismException("Heuristic model checking currently only supports simple until operators");
		}

		int bound = -1;
		
		if(exprTemp.hasBounds()) {
			bound = exprTemp.getUpperBound().evaluateInt(this.propertiesFile.getConstantValues());
			if (exprTemp.upperBoundIsStrict())
				bound--;
			if (bound < 0) {
				String boundS = exprTemp.upperBoundIsStrict() ? "<" + (bound + 1) : "<=" + bound;
				throw new PrismException("Invalid bound " + boundS + " in bounded until formula");
			}
		}
		
		mainLog.println("Starting heuristic: " + type);
		ModelExplorer modelExplorer = new CachedModelExplorer(new PrismModelExplorer(engine, modulesFile));
		//ModelExplorer modelExplorer = new PrismModelExplorer(engine, modulesFile);


		//get coalition
		Set<Integer> playersInCoalition = new HashSet<Integer>();
		if(lastExprRPATL != null) { //if this is null, we are in MDP case, I think
			int numPlayer = modulesFile.getNumPlayers();
			for(int i=0;i<numPlayer;i++) {
				Player player = modulesFile.getPlayer(i);
				if((!min && lastExprRPATL.getCoalition().contains(player.getName())) || (min && !lastExprRPATL.getCoalition().contains(player.getName()))) {
					// 0-indexed to 1-indexed
					playersInCoalition.add(i+1);
				}
			}
			//RTDP_ADJ works best and most intuitive when coalition is the maximizer; hence, if coalition is minimizing, switch it around (hence also the or in the getting of the coalition
			if(min){
				min = !min;
			}
		}



		//stateUpdate is stateUpdateSMG all the time, since I modified getStateUpdate
		StateUpdate stateUpdate = StateUpdateFactory.getStateUpdate(modulesFile.getModelType(), modelExplorer, new ThreadSafeStateValueContainer(), bound, min, epsilon);
		stateUpdate.setConstantValues(constantValues);
		stateUpdate.setCoalition(playersInCoalition);

		HeuristicNextState nextState = HeuristicNextStateFactory.getHeuristicNextState(nextStateType, modelExplorer, stateUpdate);
		Heuristic h;
		if (type == HeuristicType.RTDP_ADJ){
			h = HeuristicFactory.createHeuristic(type, this, stateUpdate, nextState, modelExplorer, min, RTDP_ADJ_optimizations, statesOA, availableOA, delta, pmin, postMax, settings);
			((HeuristicSG) h).setCoalition(playersInCoalition,modulesFile.getNumPlayers());
		}
		else {
			h = HeuristicFactory.createHeuristic(type, this, stateUpdate, nextState, modelExplorer, min);
		}
		h.setBound(bound);
		h.setVerbose(verbose);

		Expression op1 = exprTemp.getOperand1();
		if (op1 == null) {
			op1 = Expression.True();
		}
		Expression op2 = exprTemp.getOperand2();
		op1 = (Expression) op1.expandPropRefsAndLabels(propertiesFile, labelListModel);
		op1 = (Expression) op1.expandPropRefsAndLabels(propertiesFile, labelListProp);
		op2 = (Expression) op2.expandPropRefsAndLabels(propertiesFile, labelListModel);
		op2 = (Expression) op2.expandPropRefsAndLabels(propertiesFile, labelListProp);
		int operator = exprTemp.getOperator();

		Expression target = null;
		switch (operator) {
			case ExpressionTemporal.P_U://TODO: This is not how you handle until!
				throw new PrismException("Until not supported yet.");
		case ExpressionTemporal.P_F:
			target = op2;
			break;
		case ExpressionTemporal.P_G:
			target = op2;
			break;
		case ExpressionTemporal.P_W:
		case ExpressionTemporal.P_R:
		default:
			throw new PrismException("operator currently not supported for fast adaptive uniformisation");
		}
		stateUpdate.setTarget(target);

		h.compute();
		
		return new Result(new Double(h.getInitialStateValue().getLowerBound()));
	}

	private RewardStruct findRewardStruct(ExpressionReward expr) throws PrismException
	{
		RewardStruct rewStruct = null;		
		Object rs = expr.getRewardStructIndex();
		if (modulesFile == null)
			throw new PrismException("No model file to obtain reward structures");
		if (modulesFile.getNumRewardStructs() == 0)
			throw new PrismException("Model has no rewards specified");
		if (rs == null) {
			rewStruct = modulesFile.getRewardStruct(0);
		} else if (rs instanceof Expression) {
			int i = ((Expression) rs).evaluateInt(constantValues);
			rs = new Integer(i); // for better error reporting below
			rewStruct = modulesFile.getRewardStruct(i - 1);
		} else if (rs instanceof String) {
			rewStruct = modulesFile.getRewardStructByName((String) rs);
		}
		if (rewStruct == null)
			throw new PrismException("Invalid reward structure index \"" + rs + "\"");
		return rewStruct;
	}

	/**
	 * Model check an R operator.
	 */
	private Result checkExpressionReward(ExpressionReward expr) throws PrismException
	{
		mainLog.println("Starting transient probability computation using fast adaptive uniformisation...");
		PrismModelExplorer modelExplorer = new PrismModelExplorer(engine, modulesFile);
		FastAdaptiveUniformisation fau = new FastAdaptiveUniformisation(this, modelExplorer);
		ExpressionTemporal temporal = (ExpressionTemporal) expr.getExpression();
		switch (temporal.getOperator()) {
		case ExpressionTemporal.R_I:
			fau.setAnalysisType(FastAdaptiveUniformisation.AnalysisType.REW_INST);
			break;
		case ExpressionTemporal.R_C:
			fau.setAnalysisType(FastAdaptiveUniformisation.AnalysisType.REW_CUMUL);
			break;
		default:
			throw new PrismException("Currently only instantaneous or cumulative rewards are allowed.");
		}
		double time = temporal.getUpperBound().evaluateDouble(constantValues);
		RewardStruct rewStruct = findRewardStruct(expr);
		fau.setRewardStruct(rewStruct);
		fau.setConstantValues(constantValues);
		fau.computeTransientProbsAdaptive(time);
		mainLog.println("\nTotal probability lost is : " + fau.getTotalDiscreteLoss());
		mainLog.println("Maximal number of states stored during analysis : " + fau.getMaxNumStates());
		return new Result(new Double(fau.getValue()));
	}
	
}
