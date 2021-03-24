//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	* Aistis Simaitis <aistis.aimaitis@cs.ox.ac.uk> (University of Oxford)
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

package strat;

import explicit.Model;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import prism.PrismLog;

/**
 * Class to store a memoryless deterministic (MD) strategy, as a (Java) array of choice indices.
 */
public class MDStrategySparse extends MDStrategy
{
	// Model associated with the strategy
	private explicit.NondetModel model;
	// Index of choice taken in each state (wrt model above)
	// Other possible values: -1 (unknown), -2 (arbitrary), -3 (unreachable)
	private Int2IntMap choices;

	/**
	 * Creates a sparse strategy from the given choices. The choice map will be modified later on.
	 */
	public MDStrategySparse(explicit.NondetModel model, Int2IntMap choices)
	{
		this.model = model;
		this.choices = choices;
	}

	// Methods for MDStrategy

	@Override
	public int getNumStates()
	{
		return model.getNumStates();
	}

	@Override
	public boolean isChoiceDefined(int s)
	{
		return choices.get(s) >= 0;
	}

	@Override
	public void setChoiceIndex(int state, int index)
	{
		choices.put(state, index);
	}

	@Override
	public Choice getChoice(int s)
	{
		switch (choices.get(s)) {
		case -1:
			return Choice.UNKNOWN;
		case -2:
			return Choice.ARBITRARY;
		case -3:
			return Choice.UNREACHABLE;
		default:
			return Choice.INDEX;
		}
	}

	@Override
	public int getChoiceIndex(int s)
	{
		return choices.get(s);
	}

	@Override
	public Object getChoiceAction(int s)
	{
		int c = choices.get(s);
		if (c >= 0) {
			return model.getAction(s, c);
		}
		if (c == -1) {
			return "?";
		}
		if (c == -2) {
			return "*";
		}
		return "-";
	}

	// Methods for Strategy

	@Override
	public void exportInducedModel(PrismLog out)
	{
		Model dtmcInd = model.constructInducedModel(this);
		dtmcInd.exportToPrismExplicitTra(out);
	}

	@Override
	public void exportDotFile(PrismLog out)
	{
		model.exportToDotFileWithStrat(out, null, this);
	}

	@Override
	public void clear()
	{
		choices = null;
	}

	@Override public String toString()
	{
		return choices.toString();
	}
}
