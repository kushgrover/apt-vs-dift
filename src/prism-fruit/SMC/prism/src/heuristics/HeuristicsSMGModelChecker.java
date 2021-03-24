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

import explicit.FastAdaptiveUniformisation;
import explicit.ProbModelChecker.TermCrit;
import heuristics.nextstate.HeuristicNextState;
import heuristics.nextstate.HeuristicNextStateFactory;
import heuristics.nextstate.HeuristicNextStateFactory.NextState;
import heuristics.search.Heuristic;
import heuristics.search.HeuristicFactory;
import heuristics.search.HeuristicFactory.HeuristicType;
import heuristics.search.MpRtdpUnbounded;
import heuristics.update.DefaultStateValueContainer;
import heuristics.update.StateUpdate;
import heuristics.update.StateUpdateFactory;
import heuristics.update.ValueContainerFactory;
import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionPATL;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionTemporal;
import parser.ast.LabelList;
import parser.ast.ModulesFile;
import parser.ast.PropertiesFile;
import parser.ast.RewardStruct;
import prism.ModelGenerator;
import prism.ModelInfo;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismSettings;
import prism.Result;
import simulator.ModulesFileModelGenerator;
import simulator.SimulatorEngine;

/**
 * SMG model checker based on heuristic methods.
 */
public class HeuristicsSMGModelChecker extends PrismComponent
{
	// Model info
	private final ModelInfo modelInfo;
	// Model file
	private final ModulesFile modulesFile;
	// Properties file
	private final PropertiesFile propertiesFile;
	// Simulator engine
	private final SimulatorEngine engine;
	// Constants from model
	private final Values constantValues;
	// Labels from the model
	private final LabelList labelListModel;
	// Labels from the property file
	private final LabelList labelListProp;
	//Heuristic name
	private final HeuristicType type;
	//Next state type
	private final NextState nextStateType;
	//Stopping criterion
	private final double termCritParam;
	private TermCrit termCrit = TermCrit.RELATIVE;
	private final boolean verbose;
	private ExpressionTemporal currentExpression;
	private final double epsilon;

	/**
	 * Constructor.
	 */
	public HeuristicsSMGModelChecker(PrismComponent parent, ModelInfo modelInfo, ModulesFile modulesFile, PropertiesFile propertiesFile,
			SimulatorEngine engine, String heuristicName, String nextStateName, boolean verbose, double epsilon)
	{
		super(parent);
		this.modelInfo = modelInfo;
		this.modulesFile = modulesFile;
		this.propertiesFile = propertiesFile;
		this.engine = engine;

		// Get combined constant values from model/properties
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		constantValues.addValues(propertiesFile.getConstantValues());
		this.labelListModel = modulesFile.getLabelList();
		this.labelListProp = propertiesFile.getLabelList();
		this.type = HeuristicType.valueOf(heuristicName);
		this.nextStateType = NextState.valueOf(nextStateName);
		this.termCritParam = settings.getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM);
		this.verbose = verbose;
		this.epsilon = epsilon;
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
		/*
		timer = System.currentTimeMillis() - timer;
		mainLog.println("\nModel checking completed in " + (timer / 1000.0) + " secs.");

		// Print result to log
		resultString = "Result";
		if (!("Result".equals(expr.getResultName())))
			resultString += " (" + expr.getResultName().toLowerCase() + ")";
		resultString += ": " + res;
		mainLog.print("\n" + resultString + "\n");

		// Return result
		*/
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
		if (expr instanceof ExpressionPATL) {
			ExpressionPATL exprRPATL = (ExpressionPATL) expr;
			lastExprRPATL = exprRPATL;
			if (exprRPATL.getExpressionType() == ExpressionPATL.PRB) {
				return checkExpressionProb(exprRPATL.getExpressionProb());
			}
		}
		if (expr instanceof ExpressionProb) {
			ExpressionProb exprProb = (ExpressionProb) expr;
			return checkExpressionProb(exprProb);
		}
		if (expr instanceof ExpressionReward) {
			ExpressionReward exprRew = (ExpressionReward) expr;
			return checkExpressionReward(exprRew);
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
		if (expr.getRelOp().isMin()) {
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
			if (exprTemp.upperBoundIsStrict()) {
				bound--;
			}
			if (bound < 0) {
				String boundS;
				if (exprTemp.upperBoundIsStrict()) {
					boundS = "<" + (bound + 1);
				} else {
					boundS = "<=" + bound;
				}
				throw new PrismException("Invalid bound " + boundS + " in bounded until formula");
			}
		}

		mainLog.println("Starting heuristic: " + type);
		CachedModelGenerator modelGenerator = new CachedModelGenerator(new ModulesFileModelGenerator(modulesFile, this));
		StateUpdate stateUpdate = StateUpdateFactory.getStateUpdate(modulesFile.getModelType(), modelGenerator,
				ValueContainerFactory.createContainer(type, false), bound, min, epsilon);
		stateUpdate.setConstantValues(constantValues);
		if(lastExprRPATL != null) {
			throw new UnsupportedOperationException("RPATL not supported");
			// TODO: PORT
		}
		HeuristicNextState nextState = HeuristicNextStateFactory.getHeuristicNextState(nextStateType, modelGenerator, stateUpdate);
		Heuristic h = HeuristicFactory.createHeuristic(type, this, stateUpdate, nextState, modelGenerator, min);
		h.setBound(bound);
		h.setEpsilon(epsilon);
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

		final Expression target;
		switch (operator) {
		case ExpressionTemporal.P_U:
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

		// Main algorithm starts here
		h.computeProb();

		return new Result((h.getInitialStateValue().lowerBound + h.getInitialStateValue().upperBound) / 2);
	}

	private RewardStruct findRewardStruct(ExpressionReward expr) throws PrismException
	{
		RewardStruct rewardStruct = null;
		Object rewardIndex = expr.getRewardStructIndex();
		if (modulesFile == null) {
			throw new PrismException("No model file to obtain reward structures");
		}
		if (modulesFile.getNumRewardStructs() == 0) {
			throw new PrismException("Model has no rewards specified");
		}
		if (rewardIndex == null) {
			rewardStruct = modulesFile.getRewardStruct(0);
		} else if (rewardIndex instanceof Expression) {
			int i = ((Expression) rewardIndex).evaluateInt(constantValues);
			rewardIndex = i; // for better error reporting below
			rewardStruct = modulesFile.getRewardStruct(i - 1);
		} else if (rewardIndex instanceof String) {
			rewardStruct = modulesFile.getRewardStructByName((String) rewardIndex);
		}
		if (rewardStruct == null) {
			throw new PrismException("Invalid reward structure index \"" + rewardIndex + "\"");
		}
		return rewardStruct;
	}

	/**
	 * Model check an R operator.
	 */
	private Result checkExpressionReward(ExpressionReward expr) throws PrismException
	{
		if(((ExpressionTemporal) expr.getExpression()).getOperator() == ExpressionTemporal.R_M) {
			return checkExpressionMeanPayoff(expr);
		}

		mainLog.println("Starting transient probability computation using fast adaptive uniformisation...");
		ModelGenerator modelGenerator = new ModulesFileModelGenerator(modulesFile, this);
		FastAdaptiveUniformisation fau = new FastAdaptiveUniformisation(this, modelGenerator);
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
		return new Result(fau.getValue());
	}

	private Result checkExpressionMeanPayoff(ExpressionReward expr) throws PrismException
	{
		mainLog.println("Starting heuristic (mean-payoff): " + type);

		// Build rewards
		int rewardIndex = expr.getRewardStructIndexByIndexObject(modelInfo, constantValues);
		CachedModelGenerator modelGenerator = new CachedModelGenerator(new ModulesFileModelGenerator(modulesFile, this));

		int bound = -1;
		boolean min = expr.getRelOp().isMin();
		double maxRewardUpperBound = getSettings().getDouble(PrismSettings.PRISM_MDP_MP_ONTHEFLY);
		assert !Double.isNaN(maxRewardUpperBound) && 0d <= maxRewardUpperBound && maxRewardUpperBound < Double.MAX_VALUE;

		double meanPayoffPrecision = getSettings().getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM);
		double reachabilityPrecision = meanPayoffPrecision / maxRewardUpperBound;

		StateUpdate stateUpdate = StateUpdateFactory.getStateUpdate(modulesFile.getModelType(), modelGenerator,
				new DefaultStateValueContainer(), bound, min, reachabilityPrecision);
		stateUpdate.setConstantValues(constantValues);

		HeuristicNextState nextState = HeuristicNextStateFactory.getHeuristicNextState(nextStateType, modelGenerator, stateUpdate);
		MpRtdpUnbounded h = new MpRtdpUnbounded(this, stateUpdate, nextState, modelGenerator, min, meanPayoffPrecision, maxRewardUpperBound);
		h.setVerbose(verbose);
		h.setRewardIndex(rewardIndex);

		// Main algorithm starts here
		h.computeMeanPayoff();
		return new Result((h.getMeanPayoffLowerBound() + h.getMeanPayoffUpperBound()) / 2d);
	}

}
