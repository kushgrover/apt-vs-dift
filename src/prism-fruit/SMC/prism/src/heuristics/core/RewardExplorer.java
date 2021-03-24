package heuristics.core;

import explicit.MDPExplicit;
import explicit.Model;
import explicit.Product;
import explicit.rewards.MDPRewards;
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

public class RewardExplorer implements MDPRewards
{
	private final double[] stateRewards;
	private final double[][] transitionRewards;

	public RewardExplorer(ModelExplorer explorer) throws PrismException
	{
		this(explorer, 0);
	}

	public RewardExplorer(ModelExplorer explorer, int rewardIndex) throws PrismException
	{
		MDPExplicit partialModel = explorer.getPartialModel();
		ModelGenerator generator = explorer.getGenerator();

		int states = partialModel.getNumStates();
		double[] stateRewards = new double[states];
		double[][] transitionRewards = new double[states][];

		for (int stateNumber = 0; stateNumber < states; stateNumber++) {
			State state = explorer.getState(stateNumber);
			double stateReward = generator.getStateReward(rewardIndex, state);
			stateRewards[stateNumber] = stateReward;

			int choices = partialModel.getNumChoices(stateNumber);
			double[] stateTransitionRewards = new double[choices];
			for (int action = 0; action < choices; action++) {
				Object actionLabel = partialModel.getAction(stateNumber, action);
				stateTransitionRewards[action] = generator.getStateActionReward(rewardIndex, state, actionLabel);
			}
			transitionRewards[stateNumber] = stateTransitionRewards;
		}
		this.stateRewards = stateRewards;
		this.transitionRewards = transitionRewards;
	}

	@Override public double getStateReward(int stateIndex)
	{
		return stateRewards[stateIndex];
	}

	@Override public double getTransitionReward(int stateIndex, int actionIndex)
	{
		return transitionRewards[stateIndex][actionIndex];
	}

	@Override public MDPRewards liftFromModel(Product<? extends Model> product)
	{
		throw new UnsupportedOperationException("Lift not supported");
	}

	@Override public boolean hasTransitionRewards()
	{
		return true;
	}
}
