package heuristics.search.colours;

import explicit.ModelExplorer;
import parser.State;
import prism.PrismException;
import prism.PrismSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

public class ParameterDetector {

    ModelExplorer pme;
    int states=0;
    int avMax=0;
    double pmin=1d;
    int postMax=0;

    public ParameterDetector(ModelExplorer pme){
        this.pme = pme;
    }

    public void printStuff(PrismSettings settings){
        try {
            //bfs through model, remember |S|, pmin and Avmax
            LinkedList<State> worklist = new LinkedList<>();
            worklist.addLast(pme.getDefaultInitialState());

            Set<State> seen = new HashSet<>();

            while (! worklist.isEmpty()){
                State s = worklist.pop();
                seen.add(s);

                states++;
                pme.queryState(s);
                int av = pme.getNumChoices();
                avMax = av > avMax ? av : avMax;
                for(int a=0; a<av; a++){
                    int numPost = pme.getNumTransitions(a);
                    postMax = postMax > numPost ? postMax : numPost;
                    for(int t=0; t<numPost; t++){
                        double p = pme.getTransitionProbability(a,t);
                        pmin = p < pmin ? p : pmin;
                        //manage successors
                        State succ = pme.computeTransitionTarget(a,t);
                        if (!seen.contains(succ)){
                            worklist.addLast(succ);
                            seen.add(succ);
                        }
                    }
                }
                if(states%100==0){
                    System.out.println(states + ", " + seen.size() + ", " + worklist.size());
                }
            }

            System.out.println("Computation successful, here is the colourParam-String. Note that e and d are your choice.");
            System.out.println("S:"+states+";Av:"+avMax+";e:1e-3;d:0.01;p:"+pmin + ";post:"+postMax);

            try {
                String dir = settings.getString(PrismSettings.OUTPUT_DIR);
                BufferedWriter writer = new BufferedWriter(new FileWriter(dir + "/colourParams.txt"));
                writer.write("S:"+states+";Av:"+avMax+";e:1e-3;d:0.01;p:"+pmin + ";post:"+postMax);
                writer.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        catch(PrismException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

}
