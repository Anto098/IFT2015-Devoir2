package pedigree;

import java.util.ArrayList;
import java.util.Comparator;

public class PQ<Object> implements Comparator<Object> {
    private ArrayList eventHeap;
    double d;

    public PQ () {
        eventHeap = new ArrayList<Object>();
    }

    public void insert(Object object) {
        if (eventHeap.size() == 0) {
            eventHeap.add(object);
            return;
        }
        swim(object,eventHeap.size());
    }

    private void swim (Object object, int index) {
        int parentIndex = (index-1)/2;      //in java, dividing positive ints returns the whole part of the result
        while(index != 0 && compare((Object) eventHeap.get(parentIndex), object) > 0) {
            if (eventHeap.size() == index) {
                eventHeap.add(eventHeap.get(parentIndex));
            } else {
                eventHeap.set(index, eventHeap.get(parentIndex));
            }

            index = parentIndex;
            parentIndex = (index-1)/2;
        }

        if (eventHeap.size() == index) {
            eventHeap.add(object);
        } else {
            eventHeap.set(index, object);
        }
    }

    private void sink (Object v){
        int i = 0;
        int smallestChild = minChild(0);
        while (smallestChild != -1 && compare((Object) eventHeap.get(smallestChild), v) < 0 ) {
            if (eventHeap.size() == i) {
                eventHeap.add(eventHeap.get(smallestChild));
            } else {
                eventHeap.set(i, eventHeap.get(smallestChild));
            }
            i = smallestChild;
            smallestChild = minChild(i);
        }

        if (eventHeap.size() == i) {
            eventHeap.add(v);
        } else {
            eventHeap.set(i, v);
        }
    }

    private int minChild(int parentIndex) {
        int j;
        if (2*parentIndex + 2> eventHeap.size()-1 || eventHeap.get(2*parentIndex+1) == null) {      //If the parent has no children
            j=-1;
        } else if ((2*parentIndex + 2 <= eventHeap.size())
                && (compare((Object) eventHeap.get(2*parentIndex +2), (Object) eventHeap.get(2*parentIndex +1)) < 0)) {
            j=2*parentIndex+2;                          //The 2nd child (if it exists) is the smallest one
        } else {
            j=2*parentIndex+1;                          //The 1st child is the smallest one (or there is not a 2nd child)
        }
        return j;
    }

    public boolean isEmpty() {
        return eventHeap.size()==0;
    }

    public Object deleteMin() {
        if(eventHeap.size()==1){                                            // the base case if n = 1, we want to delete the only element, instead of swapping the last one (equivalent of swapping the element with itself, won't work)
            return (Object)eventHeap.remove(0);
        }
        Object min = (Object) eventHeap.get(0);
        eventHeap.set(0,eventHeap.get(eventHeap.size()-1));
        Object toSink = (Object) eventHeap.remove(eventHeap.size()-1);
        sink(toSink);
        return min;
    }

    @Override
    public int compare(Object object, Object t1) {
        if (object instanceof Sim) {
            return ((Sim) object).compareTo((Sim) t1);
        } else {
            return ((Event) object).compareTo((Event) t1);
        }
    }

    @Override
    public String toString() {
        int puissance = 0;
        for (int i=0; i<eventHeap.size();i++) {
            if (i == Math.pow(2,puissance)-1) {
                System.out.println();
                puissance++;
            }
            if(eventHeap.get(i).getClass() == Event.class){
                String type;
                Event e = (Event)eventHeap.get(i);
                if (e.type == Event.eventType.Birth) type="B";
                else if (e.type == Event.eventType.Death) type="D";
                else if (e.type == Event.eventType.Mating) type="M";
                else type="X";
                System.out.print( type+"."+e.toString()+" ");
            }
            else if (eventHeap.get(i).getClass() == Sim.class){
                System.out.print( Math.floor( ((Sim)eventHeap.get(i) ).getDeathTime())+" ");
            }

        }
        return new String();
    }

    public ArrayList getEventHeap() {
        return eventHeap;
    }

    public Sim getRandomDad(AgeModel M){
        if(eventHeap.size()==0){
            return null;
        }
        if(eventHeap.get(0).getClass()==Event.class) // if it's an Event PQ, return a new Sim (can't be a dad since it's 0 yrs old)
            return null;
        else{
            int index = (int)(M.RND.nextDouble() * eventHeap.size()-1);
            return (Sim)eventHeap.get(index);
        }
    }
}
