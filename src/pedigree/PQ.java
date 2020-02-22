package pedigree;

import java.util.ArrayList;
import java.util.Comparator;
public class PQ<Object> implements Comparator<Object> {
    private ArrayList dataHeap;
    double d;
    Types type;
    public enum Types {Events,SimDeath,SimBirth};
    public PQ (Types type) {
        this.type = type;
        dataHeap = new ArrayList<Object>();
    }

    public void insert(Object object) {
        if (dataHeap.size() == 0) {
            dataHeap.add(object);
            return;
        }
        swim(object,dataHeap.size());
    }

    private void swim (Object object, int index) {
        int parentIndex = (index-1)/2;      //in java, dividing positive ints returns the whole part of the result
        while(index > 0 && compare((Object) dataHeap.get(parentIndex), object) > 0) {
            if (dataHeap.size() == index) {
                dataHeap.add(dataHeap.get(parentIndex));
            } else {
                dataHeap.set(index, dataHeap.get(parentIndex));
            }

            index = parentIndex;
            parentIndex = (index-1)/2;
        }

        if (dataHeap.size() == index) {
            dataHeap.add(object);
        } else {
            dataHeap.set(index, object);
        }
    }

    private void sink (Object v){
        int i = 0;
        int smallestChild = minChild(0);
        while (smallestChild != -1 && compare((Object) dataHeap.get(smallestChild), v) < 0 ) {
            dataHeap.set(i, dataHeap.get(smallestChild));
            i = smallestChild;
            smallestChild = minChild(i);
        }
            dataHeap.set(i, v);
    }

    private int minChild(int parentIndex) {
        int j;
        if (2*parentIndex + 1> dataHeap.size()-1 || dataHeap.get(2*parentIndex+1) == null) {      //If the parent has no children
            j=-1;
        } else if ((2*parentIndex + 2 < dataHeap.size())
                && (compare((Object) dataHeap.get(2*parentIndex +2), (Object) dataHeap.get(2*parentIndex +1)) < 0)) {
            j=2*parentIndex+2;                          //The 2nd child (if it exists) is the smallest one
        } else {
            j=2*parentIndex+1;                          //The 1st child is the smallest one (or there is not a 2nd child)
        }
        return j;
    }

    public boolean isEmpty() {
        return dataHeap.size()==0;
    }

    public Object deleteMin() {
        if(dataHeap.size()==1){                                            // the base case if n = 1, we want to delete the only element, instead of swapping the last one (equivalent of swapping the element with itself, won't work)
            return (Object)dataHeap.remove(0);
        }
        Object min = (Object) dataHeap.get(0);
        dataHeap.set(0,dataHeap.get(dataHeap.size()-1));
        Object toSink = (Object) dataHeap.remove(dataHeap.size()-1);
        sink(toSink);
        return min;
    }

    @Override
    public int compare(Object object, Object t1) {
        if (type == Types.SimDeath) {
            return ((Sim) object).compareTo((Sim) t1);
        } else if (type == Types.Events) {
            return ((Event) object).compareTo((Event) t1);
        } else {   // if type == Types.SimBirth
            return ((Sim) object).compareBirthTime((Sim) t1);
        }
    }

    /*
    @Override
    public String toString() {
        int puissance = 0;
        for (int i=0; i<dataHeap.size();i++) {
            if (i == Math.pow(2,puissance)-1) {
                System.out.println();
                puissance++;
            }
            if(dataHeap.get(i).getClass() == Event.class){
                String type;
                Event e = (Event)dataHeap.get(i);
                if (e.type == Event.eventType.Birth) type="B";
                else if (e.type == Event.eventType.Death) type="D";
                else if (e.type == Event.eventType.Mating) type="M";
                else type="X";
                System.out.print( type+"."+e.toString()+" ");
            }
            else if (dataHeap.get(i).getClass() == Sim.class){
                System.out.print( Math.floor( ((Sim)dataHeap.get(i) ).getDeathTime())+" ");
            }

        }
        return new String();
    }
    */


    public ArrayList getDataHeap() {
        return dataHeap;
    }

    public Sim getRandomDad(AgeModel M){
        if(dataHeap.size()==0){
            return null;
        }
        if(dataHeap.get(0).getClass()==Event.class) // if it's an Event PQ, return a new Sim (can't be a dad since it's 0 yrs old)
            return null;
        else{
            int index = (int)(M.RND.nextDouble() * dataHeap.size()-1);
            return (Sim)dataHeap.get(index);
        }
    }
}
