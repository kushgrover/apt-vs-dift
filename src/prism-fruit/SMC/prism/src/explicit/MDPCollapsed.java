package explicit;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.Int2DoubleLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleSortedMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSet;
import parser.ast.ModulesFile;
import prism.Prism;
import prism.PrismException;
import prism.PrismUtils;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.function.IntConsumer;

/**
 * To store the MDP produced after the collapsing of MECs as per thesis of Luca de Alfaro
 * and the 'Reachability in MDPs: Refining Convergence of Value Iteration' paper
 *
 * @author Pranav Ashok
 */
public class MDPCollapsed extends MDPSimple
{
	int[] originalToCollapsedStateMap;

	public MDPCollapsed(MDP mdp, List<MEC> mecs, Int2DoubleMap mecRewards)
	{
		this(mdp, mecs, mecRewards, null);
	}

	/**
	 * Construct an MDPSimple object from an MDPSparse one,
	 * also collapsing the maximal end components (MECs) into a
	 * single state and applying the WR -> R transformation
	 * The last state is the '-' state and second last one is
	 * the '+' state.
	 * TODO: There are many parameters of this class  which need to be set. We are only setting trans here.
	 * TODO: Add actions. See one of the copy constructors for ref. (Reason not to add actions: we don't need rewards)
	 *
	 * @param mdp
	 * @param mecs List of MECs
	 */
	public MDPCollapsed(MDP mdp, List<MEC> mecs, Int2DoubleMap mecRewards, IntSet target)
	{
		int numOriginalStates = mdp.getNumStates();
		double maxMECReward = Collections.max(mecRewards.values());

		int numberOfCollapsedStates = 0;
		IntList transientOriginalStates = new IntArrayList();
		// For each state, map its original number to its state number in the collapsed MDP.
		// N.B. states in MECs get mapped to their representative in the collapsed MDP.
		originalToCollapsedStateMap = new int[numOriginalStates];
		Int2ObjectMap<NatBitSet> representativeMECs = new Int2ObjectAVLTreeMap<>();
		NatBitSet mecStates = NatBitSets.boundedSet(numOriginalStates);

		for (MEC mec : mecs) {
			NatBitSet statesOfMec = mec.states;
			mecStates.or(statesOfMec);
			int firstStateOfMec = statesOfMec.firstInt();
			assert firstStateOfMec != -1;
			representativeMECs.put(firstStateOfMec, statesOfMec);

			IntIterator iterator = statesOfMec.iterator();
			while (iterator.hasNext()) {
				int mecState = iterator.nextInt();
				assert originalToCollapsedStateMap[mecState] == 0;
				originalToCollapsedStateMap[mecState] = numberOfCollapsedStates;
			}

			numberOfCollapsedStates++;
		}

		IntIterator transientStateIterator = NatBitSets.complementIterator(mecStates, numOriginalStates);
		while (transientStateIterator.hasNext()) {
			int transientState = transientStateIterator.nextInt();
			transientOriginalStates.add(transientState);
			assert originalToCollapsedStateMap[transientState] == 0;
			originalToCollapsedStateMap[transientState] = numberOfCollapsedStates;
			numberOfCollapsedStates++;
		}

		// All transient + all collapsed MECs + plus + minus
		super.numStates = numberOfCollapsedStates + 2;

		// The first idx numStates are going to the union of the transient states and the representatives of MECs
		// The (numStates-1)th state is going to be the '+' state and (numStates)th state is going to be '-'
		super.trans.clear();
		for (int i = 0; i < numStates; i++) {
			// Will be replaced later
			super.trans.add(new Distribution[] {});
		}

		// Go through each original state and ...
		IntIterator transientOriginalStateIterator = transientOriginalStates.iterator();
		while (transientOriginalStateIterator.hasNext()) {
			int transientState = transientOriginalStateIterator.nextInt();
			// This state is a transient one, map its distribution
			// Get every distributions out of it and add it to trans
			int stateChoices = mdp.getNumChoices(transientState);
			List<Distribution> distributions = new ArrayList<>(stateChoices);
			for (int choiceNumber = 0; choiceNumber < stateChoices; choiceNumber++) {
				Distribution distribution = new Distribution();
				for (Int2DoubleMap.Entry transitionEntry : mdp.getTransitions(transientState, choiceNumber)) {
					// Destination in the original MDP
					// Destination in the collapsed MDP
					int collapsedDestination = originalToCollapsedStateMap[transitionEntry.getIntKey()];
					distribution.add(collapsedDestination, transitionEntry.getDoubleValue());
				}
				distributions.add(distribution);
			}
			int collapsedTransientState = originalToCollapsedStateMap[transientState];
			super.trans.set(collapsedTransientState, distributions.toArray(new Distribution[distributions.size()]));
		}

		for (Int2ObjectMap.Entry<NatBitSet> mecEntry : representativeMECs.int2ObjectEntrySet()) {
			// This state is a representative of a MEC
			// Find all distributions which leave the MEC and give them to s
			int mecRepresentative = mecEntry.getIntKey();
			int mecCollapsedState = originalToCollapsedStateMap[mecRepresentative];
			NatBitSet mec = mecEntry.getValue();
			List<Distribution> choices = new ArrayList<>();

			// TODO Use the values from MEC class here?
			// Iterate over all states in the MEC
			IntIterator mecStateIterator = mec.iterator();
			while (mecStateIterator.hasNext()) {
				int j = mecStateIterator.nextInt();
				// And their choices
				for (int i = 0; i < mdp.getNumChoices(j); i++) {
					if (mdp.allSuccessorsInSet(j, i, mec)) {
						continue;
					}
					// If state j has any action i which has an edge going outside the MEC, compute this action/distribution
					Distribution distribution = new Distribution();
					double outgoingProbability = 0d;
					for (Int2DoubleMap.Entry entry : mdp.getTransitions(j, i)) {
						double prob = entry.getDoubleValue();
						int collapsedDestination = originalToCollapsedStateMap[entry.getIntKey()];
						if (mecCollapsedState == collapsedDestination) {
							// Self-loop in the collapsed MEC can be omitted for reachability
							continue;
						}
						// If the distribution does not contain representative, 0.0 is returned
						distribution.add(collapsedDestination, distribution.get(collapsedDestination) + prob);
						outgoingProbability += prob;
					}
					Distribution Distribution;
					if (outgoingProbability < 1d) {
						// Rescale the distribution
						double OutgoingProbability = outgoingProbability;
						Distribution = new Distribution();
						distribution.forEach((key, value) -> Distribution.add(key, value / OutgoingProbability));
					} else {
						Distribution = distribution;
					}
					choices.add(Distribution);
				}
			}

			// Since s represents a MEC, apply the reachability transformation here
			Distribution toPlusAndMinus = new Distribution();
			assert mecRewards.containsKey(mecRepresentative);
			double reward = mecRewards.get(mecRepresentative) / maxMECReward;

			// Transition to + state with g(s)/gmax probability
			if (reward > 0.) {
				toPlusAndMinus.add(numStates - 2, reward);
			}

			// Transition to + state with 1 - g(s)/gmax probability
			if (reward < 1.) {
				toPlusAndMinus.add(numStates - 1, 1d - reward);
			}

			choices.add(toPlusAndMinus);
			super.trans.set(mecCollapsedState, choices.toArray(new Distribution[choices.size()]));
		}

		// In the case of interval iteration
		if (target != null) {
			// For each target state, add prob 1 transition to + state
			target.forEach((IntConsumer) targetState -> {
				int collapsedTargetState = originalToCollapsedStateMap[targetState];
				super.trans.set(collapsedTargetState, new Distribution[] { new Distribution(numStates - 2, 1d) });
			});
			// If we don't have this, there could be a case when a transient state is a
			// a goal state - but it doesn't have a transition to +.
		}

		// Add self loop to + state
		Distribution plusLoopDistribution = new Distribution(numStates - 2, 1d);
		super.trans.set(numStates - 2, new Distribution[] { plusLoopDistribution });

		// Add self loop to - state
		Distribution minusLoopDistribution = new Distribution(numStates - 1, 1d);
		super.trans.set(numStates - 1, new Distribution[] { minusLoopDistribution });

		mdp.getInitialStates().forEach((IntConsumer) initialState -> {
			int collapsedInitialState = originalToCollapsedStateMap[initialState];

			if (!this.initialStates.contains(collapsedInitialState)) {
				this.initialStates.add(collapsedInitialState);
			}
		});

		assert (trans.size() == numStates);
		assert super.trans.stream().allMatch(Objects::nonNull);
	}

	public int getCollapsedStateIndexForOriginalState(int originalState)
	{
		return originalToCollapsedStateMap[originalState];
	}

	public int getGoalState()
	{
		return getNumStates() - 2;
	}

	public int getTrapState()
	{
		return getNumStates() - 1;
	}

	/**
	 * Converts a collapsed MDP into a PRISM language file string.
	 *
	 * @return The modules file as a String
	 */
	private String getAsPrismLanguage()
	{
		StringBuilder outputBuilder = new StringBuilder(1000);
		outputBuilder.append(getModelType().keyword()).append('\n');
		outputBuilder.append("module M\nx : [0..").append(numStates - 1).append("] init ").append(initialStates.getInt(0)).append(";\n");

		Int2DoubleSortedMap sorted = new Int2DoubleLinkedOpenHashMap();
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < getNumChoices(i); j++) {
				// Extract transitions and sort by destination state index (to match PRISM-exported files)
				for (Int2DoubleMap.Entry entry : getTransitions(i, j)) {
					if (entry.getDoubleValue() != 0.0d) {
						sorted.put(entry.getIntKey(), entry.getDoubleValue());
					}
				}

				// Print out (sorted) transitions
				Object action = getAction(i, j);
				if (action == null) {
					outputBuilder.append("[]");
				} else {
					outputBuilder.append('[').append(action).append(']');
				}
				outputBuilder.append("x=").append(i).append("->");

				boolean first = true;
				for (Int2DoubleMap.Entry e : sorted.int2DoubleEntrySet()) {
					if (first) {
						first = false;
					} else {
						outputBuilder.append("+");
					}
					// Note use of PrismUtils.formatDouble to match PRISM-exported files
					outputBuilder.append(PrismUtils.formatDouble(e.getDoubleValue())).append(":(x'=").append(e.getIntKey()).append(')');
				}
				outputBuilder.append(";\n");
				sorted.clear();
			}
		}
		outputBuilder.append("endmodule\n");
		return outputBuilder.toString();
	}

	public ModulesFile generateModulesFile() throws PrismException
	{
		String mdpAsString = this.getAsPrismLanguage();
		try (InputStream inputStream = new ByteArrayInputStream(mdpAsString.getBytes(StandardCharsets.UTF_8))) {
			ModulesFile mf = Prism.getPrismParser().parseModulesFile(inputStream);
			Prism.releasePrismParser();
			return mf;
		} catch (InterruptedException e) {
			throw new PrismException("Error with fetching parser/parser in use");
		} catch (IOException e) {
			throw new PrismException("Can't happen");
		}
	}

	@Override public Set<Info> getMDPInformation()
	{
		return EnumSet.of(Info.NO_TRANSIENT_MEC);
	}
}
