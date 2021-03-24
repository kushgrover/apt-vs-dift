package heuristics.core;

import explicit.MDPExplicit;
import explicit.MDPModelChecker;
import explicit.ModelCheckerResult;
import explicit.StateValues;
import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionTemporal;
import parser.ast.ModulesFile;
import parser.ast.PropertiesFile;
import prism.ModelInfo;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismSettings;
import prism.Result;
import simulator.ModulesFileModelGenerator;

import java.util.concurrent.TimeUnit;

import static com.google.common.base.Preconditions.checkArgument;

public class CoreModelChecker extends PrismComponent
{
	private final ModelInfo modelInfo;
	private final double epsilon;
	private final Values constantValues;
	private final ModulesFileModelGenerator model;
	private InfiniteCoreLearner infiniteCoreLearner = null;
	private FiniteCoreLearner finiteCoreLearner = null;

	public CoreModelChecker(PrismComponent parent, ModelInfo modelInfo, ModulesFile modulesFile, PropertiesFile propertiesFile, double epsilon)
			throws PrismException
	{
		super(parent);
		this.modelInfo = modelInfo;
		this.epsilon = epsilon;
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		constantValues.addValues(propertiesFile.getConstantValues());
		model = new ModulesFileModelGenerator(modulesFile, this);
	}

	private MDPExplicit computeInfiniteEpsilonCore(double requiredPrecision) throws PrismException
	{
		if (infiniteCoreLearner == null) {
			infiniteCoreLearner = new InfiniteCoreLearner(this, model);
		}
		return infiniteCoreLearner.learn(requiredPrecision);
	}

	private MDPExplicit computeStepBoundedEpsilonCore(int stepBound, double requiredPrecision) throws PrismException
	{
		if (finiteCoreLearner == null) {
			finiteCoreLearner = new FiniteCoreLearner(this, model);
		}
		return finiteCoreLearner.learn(stepBound, requiredPrecision);
	}

	/**
	 * Model check a property.
	 */
	public Result check(Expression expr) throws PrismException
	{
		// TODO Add more expressions like reachability and step bounded properties
		checkArgument(expr instanceof ExpressionReward, "Only reward properties are supported");
		ExpressionReward rewardExpr = (ExpressionReward) expr;

		checkArgument(((ExpressionTemporal) rewardExpr.getExpression()).getOperator() == ExpressionTemporal.R_M, "Only mean payoff property is supported");
		long timer = System.nanoTime();

		Result res = checkExpressionMeanPayoff(rewardExpr);

		// Model checking complete
		timer = System.nanoTime() - timer;
		mainLog.println(String.format("%nModel checking completed in %f secs.", timer / (double) TimeUnit.SECONDS.toNanos(1)));

		return res;
	}

	private Result checkExpressionMeanPayoff(ExpressionReward expr) throws PrismException
	{
		int rewardIndex = expr.getRewardStructIndexByIndexObject(modelInfo, constantValues);

		boolean min = expr.getRelOp().isMin();
		double maxRewardUpperBound = getSettings().getDouble(PrismSettings.PRISM_MDP_MP_ONTHEFLY);
		assert !Double.isNaN(maxRewardUpperBound) && 0d <= maxRewardUpperBound && maxRewardUpperBound < Double.MAX_VALUE;
		double reachabilityPrecision = epsilon / maxRewardUpperBound;

		MDPExplicit epsilonCore = computeInfiniteEpsilonCore(reachabilityPrecision);
		// TODO Don't use the MDPMPChecker directly here, adapt the code to reuse the gained knowledge
		MDPModelChecker mdpmc = new MDPModelChecker(this);
		RewardExplorer rewards = new RewardExplorer(infiniteCoreLearner.getCollapsingExplorer().getPartialExplorer(), rewardIndex);

		ModelCheckerResult result =  mdpmc.computeAverageReward(epsilonCore, rewards, min);
		Result res = new Result();
		res.setVector(StateValues.createFromDoubleArray(result.soln, epsilonCore));
		return res;
	}
}
