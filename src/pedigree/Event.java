package pedigree;

public class Event implements Comparable<Event> {
    double time;
    Sim subject;
    eventType type;
    public enum eventType {Birth, Death, Mating};

    public Event (double time, Sim subject,eventType type) {
        this.time = time;
        this.subject = subject;
        this.type = type;
    }

    public int compareTo(Event event) {
        return (time >= event.time) ? 1 : -1;   //If the 2 events occur at the same time, then we don't really care and we return that "this" event is the first one
    }

    @Override
    public String toString() {
        return ""+(int)Math.floor(time);
    }
}
