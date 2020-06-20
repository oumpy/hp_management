Title: Contact
Date: 2020.03.13
Modified: 2020.03.17
Summary: Python会への連絡用メールフォーム

<!-- ![]({attach}images/computer-1209641_960_7201.jpg) -->

Python会に興味のある方、解析依頼、コンペ出場依頼などありましたら、E-mail, Twitter等にてお気軽にご連絡ください。

<form id="fs-frm" name="simple-contact-form" accept-charset="utf-8" action="https://formspree.io/xzbvrlev" method="post">
  <fieldset id="fs-frm-inputs">
    <label for="full-name">Name (必須)</label><br>
    <input type="text" name="name" id="full-name" placeholder="" required=""><br><br>
    <label for="email-address">Email Address (必須)</label><br>
    <input type="email" name="replyto" id="email-address" placeholder="" required=""><br><br>
    <label for="message">Message (必須)</label><br>
    <textarea rows="5" cols="100" name="message" id="message" placeholder="" required=""></textarea>
    <input type="hidden" name="subject" id="email-subject" value="Contact Form Submission">
  </fieldset>
  <input type="submit" value="Submit">
</form>

<style>
input{
	background-color:#f5f5f5;
	width:400px;
	max-width: 95%;
}

textarea{
	background-color:#f5f5f5;
	max-width:95%;
}
</style>
