Set2b=brewer.pal(n = 8, name = "YlGnBu")

hist(d_cap$ageatcap)

names(d_cap)
d_cap$Sex <- as.factor(d_cap$sex)
d_cap$Sex <- factor(d_cap$sex,labels=c("Male","Female"))
ageatcap_plot <- ggplot(d_cap) + geom_histogram(aes(x=ageatcap,fill=Sex))+
                facet_wrap(.~Sex)+
                scale_fill_manual(values=Set2b[c(5,3)])+
                theme_bw()+
                ylab("Count")+
                xlab("Age at Capture (Years)")+
                ggtitle("Age at Capture of GPS Collared Deer")+
                theme(axis.text=element_text(size=12),
                axis.title=element_text(size=14,face="bold"),
                title =element_text(size=12, face='bold'),
                strip.text.x = element_text(size = 12))
ageatcap_plot
ggsave("figures/ageatcap_plot.png",ageatcap_plot,height=4,width=7)
