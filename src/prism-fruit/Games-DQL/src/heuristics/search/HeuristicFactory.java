//==============================================================================
//	
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Mateusz Ujma <mateusz.ujma@cs.ox.ac.uk> (University of Oxford)
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

package heuristics.search;

import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.nextstate.HeuristicNextStateFactory;
import heuristics.search.colours.Black;
import heuristics.search.colours.Grey;
import heuristics.search.colours.ParameterDetector;
import heuristics.search.colours.White;
import heuristics.update.StateUpdate;
import explicit.ModelExplorer;
import prism.PrismException;
import prism.PrismSettings;

public class HeuristicFactory
{
	public enum HeuristicType {
		RTDP, LRTDP, RTDP_PARALLEL, LRTDP_PARALLEL, ILAO, RTDP_UNBOUNDED, RTDP_ADJ
	}
	
	public static Heuristic createHeuristic(HeuristicType type, HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException{
		switch(type) {
			case RTDP:
				return new HeuristicRTDP(mc, su, nextState, pme, min);
			case LRTDP:
				return new HeuristicLRTDP(mc, su, nextState, pme, min);
			case RTDP_PARALLEL:
				return new HeuristicRTDPParallel(mc, su, nextState, pme, min);
			case LRTDP_PARALLEL:
				return new HeuristicLRTDPParallel(mc, su, nextState, pme, min);
			case ILAO:
				return new HeuristicILAO(mc, su, nextState, pme, min);
			//case RTDP_UNBOUNDED:
			//	return new HeuristicRTDP_Unbounded(mc, su, nextState, pme, min);
			case RTDP_ADJ: //All the new heuristics, doing collapsing and deflating and black/grey/white setting
				throw new PrismException("Not possible in this version of createHeuristic; need to give opts to the method and use different call.");

		}
		throw new PrismException("Unimplemented heuristic: " + type);
	}

	public static Heuristic createHeuristic(HeuristicType type, HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min, int opts, int statesOA, int availableOA, double delta, double pmin, int postMax, PrismSettings settings) throws PrismException{
		if (opts % 8 != 0){
			if(opts==3){
				//special case of my parameter detector
				ParameterDetector pd = new ParameterDetector(pme);
				pd.printStuff(settings);
				System.exit(0);
			}
			//not white setting, so nextState has to be HIGH_PROB to allow sampling
			nextState = HeuristicNextStateFactory.getHeuristicNextState(HeuristicNextStateFactory.NextState.HIGH_PROB, pme, su);
			try{
				double test = statesOA + availableOA + delta + pmin;
			}
			catch (Exception e){
				throw new PrismException("Grey or Black requested, but colourParams not set.");
			}
		}
		//Note: stateUpdate is stateUpdateSMG; since in (Strg+click on createHeuristic) I modified getStateUpdate to do that
		HeuristicSG h;
		switch(type) {
			case RTDP_ADJ: //All the new heuristics, doing collapsing and deflating and black/grey/white setting
				switch(opts%8){
					//cases for different settings/versions
					case 0:
						h = new White(mc,su,nextState,pme,min);
						break;
					case 1:
						h = new Grey(mc,su,nextState,pme,min, statesOA, availableOA, delta,pmin, postMax);
						break;
					case 2:
						h = new Black(mc,su,nextState,pme,min, statesOA, availableOA, delta,pmin, postMax);
						break;
//					case 3:
//						h = new Old(mc,su,nextState,pme,min);
//						break;
//					case 4:
//						h = new DQL_compareBstat(mc,su,nextState,pme,min);
//						break;
//					case 5:
//						h = new DQL_compareBold(mc,su,nextState,pme,min);
//						break;
//					case 6:
//						h = new DQL_compareU(mc,su,nextState,pme,min);
//						break;
//					case 7:
					default:
						throw new PrismException("Unknown version of HeuristicSG (version=" + (opts%8) + ")");
				}
				h.setOpts(opts);
				return h;
			default:
				return createHeuristic(type,mc,su,nextState,pme,min);

		}
	}
}
