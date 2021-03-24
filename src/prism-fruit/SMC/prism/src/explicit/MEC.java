/*
 * Copyright (C) 2016  (Salomon Sickert)
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

package explicit;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIterator;

import java.util.Objects;

/**
 * MEC describes a MEC for a given model. It does not store the model itself. Only the states belonging to the MEC and
 * the allowed actions.
 */
public final class MEC
{
	public final NatBitSet states;
	public final Int2ObjectMap<NatBitSet> actions;

	private MEC(NatBitSet s, Int2ObjectMap<NatBitSet> a)
	{
		states = s;
		actions = a;
	}

	public static MEC createMEC(MDPSimple model, NatBitSet states)
	{
		Int2ObjectMap<NatBitSet> actions = new Int2ObjectOpenHashMap<>(states.size());
		boolean changed = true;

		while (changed) {
			changed = false;

			IntIterator stateIterator = states.iterator();
			while (stateIterator.hasNext()) {
				int i = stateIterator.nextInt();
				NatBitSet act = actions.get(i);

				Distribution[] distributions = model.trans.get(i);

				if (distributions != null) {
					if (act == null) {
						act = NatBitSets.boundedFilledSet(distributions.length);
						actions.put(i, act);
					}

					IntIterator actIterator = act.iterator();
					while (actIterator.hasNext()) {
						int j = actIterator.nextInt();
						if (!distributions[j].isSubsetOf(states)) {
							act.clear(j);
						}
					}
				}

				if (act == null || act.isEmpty()) {
					stateIterator.remove();
					actions.remove(i);
					changed = true;
				}
			}
		}

		return new MEC(states, actions);
	}

	public static MEC createMEC(NondetModel model, NatBitSet states)
	{
		if (model instanceof MDPSimple) {
			return createMEC((MDPSimple) model, states);
		}

		Int2ObjectMap<NatBitSet> actions = new Int2ObjectOpenHashMap<>(states.size());
		boolean changed = true;

		while (changed) {
			changed = false;

			IntIterator stateIterator = states.iterator();
			while (stateIterator.hasNext()) {
				int state = stateIterator.nextInt();

				NatBitSet stateActions = actions.get(state);
				if (stateActions == null) {
					int numChoices = model.getNumChoices(state);
					stateActions = NatBitSets.boundedFilledSet(numChoices, numChoices);
					actions.put(state, stateActions);
				}

				IntIterator actionIterator = stateActions.iterator();
				while (actionIterator.hasNext()) {
					int action = actionIterator.nextInt();
					if (!model.allSuccessorsInSet(state, action, states)) {
						stateActions.clear(action);
					}
				}

				if (stateActions.isEmpty()) {
					states.clear(state);
					actions.remove(state);
					changed = true;
				}
			}
		}

		return new MEC(states, actions);
	}

	@Override
	public boolean equals(Object o)
	{
		if (this == o)
			return true;
		if (o == null || getClass() != o.getClass())
			return false;
		MEC lightMEC = (MEC) o;
		return Objects.equals(states, lightMEC.states);
	}

	@Override public String toString()
	{
		return states.toString();
	}

	@Override
	public int hashCode()
	{
		return Objects.hash(states);
	}

	public int size()
	{
		return states.size();
	}
}
