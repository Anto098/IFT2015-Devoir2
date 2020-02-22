/*
 * Copyright 2020 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package pedigree;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;

/**
 *
 * Gompertz-Makeham distribution for lifespan.
 *
 * Parametrized by accident, death rate and scale.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class AgeModel
{
    private double stable_rate;
    Random RND = new Random();
    private final double death_rate;
    private final double accident_rate;
    private final double age_factor;

    private static final double DEFAULT_ACCIDENT_RATE = 0.01; // 1% chance of dying per year
    private static final double DEFAULT_DEATH_RATE = 12.5;
    private static final double DEFAULT_SCALE = 100.0; // "maximum" age [with death rate 1]

    public AgeModel(double accident_rate, double death_rate, double age_scale) {
        this.death_rate = death_rate;
        this.age_factor = Math.exp(age_scale/death_rate);
        this.accident_rate = accident_rate;

    }

    /**
     * Instantiation with default values (human).
     */
    public AgeModel() {
        this(DEFAULT_ACCIDENT_RATE, DEFAULT_DEATH_RATE, DEFAULT_SCALE);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(getClass().getName());
        sb.append("[acc ").append(accident_rate).append(", age ").append(death_rate).append(", agefactor ").append(age_factor);
        sb.append("]");
        return sb.toString();
    }

    /**
     * Probability of surviving past the given age
     *
     * @param age
     * @return probability of dying after the given age
     */
    public double getSurvival(double age) {
        return Math.exp(-accident_rate*age -death_rate*Math.expm1(age/death_rate)/age_factor);
    }

    /**
     * Expected time span (TS) for mating: average number of children will be TS/matingrate.
     *
     * @param min_age minimum age of sexual maturity
     * @param max_age maximum age of parenting
     * @return
     */
    public double expectedParenthoodSpan(double min_age, double max_age) {

        // integration of the survival function over the mating age

        // numerical integration with simpson's rule, dynamic setting of resolution

        int n = 1; // number of intervals along the range
        double d = (max_age-min_age)/n;

        double st = d*0.5*(getSurvival(min_age)+getSurvival(max_age));


        double espan = 0.0;
        double old_espan = -1.0; // does not matter much

        for (int iter=1; iter<20;iter++)
        {
            double x0=min_age+d*0.5;
            double s2=0.0;
            for (int i=0;i<n;i++)
            {
                double x = x0+i*d;
                s2 += getSurvival(x);
            }
            double old_st = st;
            st = 0.5*(st+d*s2); // simple trapezoidal
            espan = (4.0*st-old_st)/3.0; // Simpson's ... better than st

            n = n*2;
            d=d*0.5;

            if (iter>5 // first five iteration kept
                    && (Math.abs(old_espan-espan)<1e-7*old_espan
                    || (espan==0.0 && old_espan==0.0) ))
                break;
            old_espan = espan;
        }
        return espan;
    }

    /**
     * Exponentially distributed random variable.
     *
     *
     * @param RND random number generator
     * @param rate inverse of the mean
     * @return Exponential(rate)
     */

    public static double randomWaitingTime(Random RND, double rate)
    {
        return -Math.log(RND.nextDouble())/rate;
    }

    /**
     * A random value with the specified lifespan distribution.
     *
     * @param RND Pseudorandom number generator for uniform[0,1]
     *
     * @return a random value distributed by Gomperz-Makeham
     */
    public double randomAge(Random RND) {
        // pseudorandom by exponential for accident-related death
        double accidental_death = -Math.log(RND.nextDouble())/accident_rate;
        // pseudorandom by Gompertz for old-age
        double u = RND.nextDouble();
        double age_death = death_rate*Math.log1p(-Math.log(u)/death_rate*age_factor);

        return Math.min(age_death, accidental_death);
    }


    /**
     * To simulate from an initial founders population
     * @param n         size of the founders population
     * @param Tmax      maximal time for the simulation
     */
    double Time = 0;
    PQ<Event> eventQ;
    PQ<Sim> simsQ;
    double totalWaiting =0;
    int totalMatings =0;
    private void simulate(int n, double Tmax) {
        int lastCentury=0;
        System.out.println(0+","+n);

        eventQ = new PQ<>(PQ.Types.Events);  // file de priorité
        simsQ = new PQ<>(PQ.Types.SimDeath);      // simsQ
        for (int i = 0; i < n; i++) {
            Sim fondateur = new Sim(Sim.getRandomSex()); // sexe au hasard, naissance à 0.0
            Event E = new Event(0,fondateur,Event.eventType.Birth);
            eventQ.insert(E); // insertion dans la file de priorité
        }

        while (!eventQ.isEmpty()) {

            if (Math.floor(Time/100) > lastCentury) {
                lastCentury++;
                System.out.println((int)Time+","+simsQ.getDataHeap().size());
            }

            Event E = eventQ.deleteMin(); // prochain événement
            Time = E.time;                // update time
            if (E.time > Tmax) break;     // arrêter à Tmax
            if (E.subject.getDeathTime() >= E.time){
                if(E.type == Event.eventType.Birth){                         // BIRTH
                    if (!E.subject.isFounder()) E.subject.getMother().nbBabies++;

                    double deathTime = Time+randomAge(RND);                                                 // make death time
                    eventQ.insert(new Event(deathTime,E.subject,Event.eventType.Death));                    // if it's a birth, we add a deathEvent at Time + age of the sim
                    E.subject.setDeathTime(deathTime);                                                      // set death time
                    simsQ.insert(E.subject);                                                                // insert the sim in the simsQ
                    if(E.subject.getSex() == Sim.Sex.F){
                        double waitingTime = randomWaitingTime(RND,stable_rate);
                        totalWaiting+=waitingTime;
                        totalMatings++;
                        eventQ.insert(new Event(Time+waitingTime,E.subject,Event.eventType.Mating)); // if subject is girl, create new mating event
                    }
                    if(!E.subject.isFounder() && E.subject.getMother().isMatingAge(Time)){                  // After giving birth, can the mom still have another child?
                        double waitingTime = randomWaitingTime(RND,stable_rate);
                        double nextMatingTime = Time+waitingTime;
                        if(nextMatingTime<E.subject.getMother().getDeathTime()){                            // if mother will mate again before she dies, we add the event
                            totalWaiting+=waitingTime;
                            totalMatings++;
                            eventQ.insert(new Event(nextMatingTime,E.subject.getMother(),Event.eventType.Mating));
                        }                                                                                   // else we don't add the event because it will only increase the size of our eventQ but will never happen
                    }
                }
                else if(E.type == Event.eventType.Death) {                   // DEATH
                    Sim deadSim = simsQ.deleteMin();                                                                      // remove the sim from the list of the population
                }
                else if(E.type == Event.eventType.Mating) {                  // MATING

                    double age = Time - E.subject.getBirthTime();
                    if (Sim.MIN_MATING_AGE_F <= age && age <= Sim.MAX_MATING_AGE_F) {                       // if subject is between min age and max age

                        Sim papa = E.subject.chooseMate(this, simsQ);                                       // choose new dad
                        if(papa == null){
                            double waitingTime = randomWaitingTime(RND,stable_rate);
                            totalWaiting+=waitingTime;
                            totalMatings++;
                            eventQ.insert(new Event(Time+waitingTime,E.subject,Event.eventType.Mating)); // if she tried too many times (it's too hard to find a mate for her right now), we make a new mating event for her

                            continue;
                        }
                        Sim baby = new Sim(E.subject, papa, Time, Sim.getRandomSex());         // make new baby
                        Event birth = new Event(Time, baby, Event.eventType.Birth);            // make new birth event
                        eventQ.insert(birth);                                                               // put birth event in eventQ

                        E.subject.setMate(papa);
                        papa.setMate(E.subject);
                    } else if (age <= Sim.MIN_MATING_AGE_F) {                                               // ADDED : si fille trop jeune, remettre un mating plus tard
                        double waitingTime = randomWaitingTime(RND,stable_rate);
                        totalWaiting+=waitingTime;
                        totalMatings++;
                        eventQ.insert(new Event(Time + waitingTime,E.subject,Event.eventType.Mating));
                    }
                }
            }

        }
        System.out.println((int)Time+","+simsQ.getDataHeap().size());
        getCoalescence(simsQ, Tmax);
    }

    private void getCoalescence(PQ<Sim> simsQ, double Tmax){
        PQ<Sim> hCoaQ = new PQ<>(PQ.Types.SimBirth);        // new Min Heap to put all our alive people in
        PQ<Sim> fCoaQ = new PQ<>(PQ.Types.SimBirth);        // H F = Homme Femme

        HashSet<Integer> fAncestors = new HashSet<>();              // to keep all the ancestors that we've encountered
        HashSet<Integer> hAncestors = new HashSet<>();

        ArrayList<Sim> simsToAdd = new ArrayList<>();

        while(!simsQ.isEmpty()) {                           //Preliminary check to see if some parents are still alive. If they are, we don't add them now to the CoaQ because it would duplicate them later
            Sim simToAdd = simsQ.deleteMin();
            Sim parent = (simToAdd.getSex() == Sim.Sex.F) ? simToAdd.getMother() : simToAdd.getFather();
            if(simToAdd.getSex() == Sim.Sex.F) {            // if it's a girl
                simsToAdd.add(simToAdd);
                if (parent.getDeathTime() > Tmax) {         //if her parent is still alive, we add the parent to the ancestors
                    fAncestors.add(parent.getSimIdent());
                }
            } else {                                        // if it's a boy
                simsToAdd.add(simToAdd);
                if (parent.getDeathTime() > Tmax) {
                    hAncestors.add(parent.getSimIdent());   //if his parent is still alive, we add the parent to the ancestors
                }
            }
        }

        for (Sim simToAdd : simsToAdd) {                    // Put sims in the HCoaQ and FCoaQ
            if(simToAdd.getSex() == Sim.Sex.F) {            // if it's a girl
                if (!fAncestors.contains(simToAdd.getSimIdent())) {
                    fCoaQ.insert(simToAdd);
                }
            } else {                                        // if it's a boy
                if (!hAncestors.contains(simToAdd.getSimIdent())) {
                    hCoaQ.insert(simToAdd);
                }
            }
        }

        fAncestors = new HashSet<>();              // we reset the ancestor sets because we used it for our preliminary check
        hAncestors = new HashSet<>();

        HashMap<Double,Integer> hCoalPoints = new HashMap<>();
        HashMap<Double,Integer> fCoalPoints = new HashMap<>();  // keep all the coalescence points

        System.out.println("\n\n\n\n\n\n\n");
        treatSimCoalescence(hCoaQ,hAncestors,hCoalPoints);
        System.out.println("\n\n\n\n\n\n\n");
        treatSimCoalescence(fCoaQ,fAncestors,fCoalPoints);
    }

    private void treatSimCoalescence(PQ<Sim> coaQ, HashSet<Integer> ancestors, HashMap<Double,Integer> coalPoints) {
        while(!coaQ.isEmpty()) {
            Sim sim = coaQ.deleteMin();

            if (!sim.isFounder()) {
                Sim parent = (sim.getSex() == Sim.Sex.F) ? sim.getMother() : sim.getFather();
                if (ancestors.contains(parent.getSimIdent())) {                             // if the parent is already in the HashSet of ancestors
                    coalPoints.put(sim.getBirthTime(),coaQ.getDataHeap().size());                // we add a coalpoint with his child's birthtime
                    System.out.println(""+sim.getBirthTime()+","+coaQ.getDataHeap().size());
                } else {                                                                    // else the parent is not already in the HashSet of ancestors
                    ancestors.add(parent.getSimIdent());                                          // add him to ancestors
                    coaQ.insert(parent);                                                          // add him to coaQ
                }
            }
        }
    }

    /**
     * Test for tabulating random lifespans from command line.
     *
     * @param args accident-rate death-rate [scale]
     */
    public static void main(String[] args) {
        int arg_idx = 0;
        AgeModel M;
        if (args.length > 1) {
            int n = Integer.parseInt(args[arg_idx++]);
            double Tmax = Double.parseDouble(args[arg_idx]);
            M = new AgeModel();
            double span = M.expectedParenthoodSpan(Sim.MIN_MATING_AGE_F, Sim.MAX_MATING_AGE_F);
            M.stable_rate = 2/span;
            M.simulate(n,Tmax);

        } else {
            System.out.println("Please enter arguments in the form {n, Tmax}");
        }

    }

}
